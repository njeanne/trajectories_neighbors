#!/usr/bin/env python3

"""
Created on 12 Feb. 2025
"""

__author__ = "Nicolas JEANNE"
__copyright__ = "GNU General Public License"
__email__ = "jeanne.n@chu-toulouse.fr"
__version__ = "1.0.0"

import argparse
import logging
import os
import pickle
import re
import sys
import yaml

from mpi4py import MPI
import pandas as pd
import pytraj as pt


def restricted_float(value_to_inspect):
    """Inspect if a float is between 0.0 and 100.0

    :param value_to_inspect: the value to inspect
    :type value_to_inspect: str
    :raises ArgumentTypeError: is not between 0.0 and 100.0
    :return: the float value if float_to_inspect is between 0.0 and 100.0
    :rtype: float
    """
    x = float(value_to_inspect)
    if x < 0.0 or x > 100.0:
        raise argparse.ArgumentTypeError(f"{x} not in range [0.0, 100.0]")
    return x


def restricted_positive(value_to_inspect):
    """Inspect if the value is positive.

    :param value_to_inspect: the value to inspect
    :type value_to_inspect: str
    :raises ArgumentTypeError: is not > 0.0
    :return: the float value if float_to_inspect is between 0.0 and 100.0
    :rtype: float
    """
    x = float(value_to_inspect)
    if x < 0.0:
        raise argparse.ArgumentTypeError(f"{x} not a positive value.")
    return x


def restricted_angle(value_to_inspect):
    """Inspect if an angle value is between 0 and 359.

    :param value_to_inspect: the value to inspect
    :type value_to_inspect: str
    :raises ArgumentTypeError: is not between 0.0 and 100.0
    :return: the float value if float_to_inspect is between 0.0 and 100.0
    :rtype: float
    """
    x = int(value_to_inspect)
    if x < 0 or x > 359:
        raise argparse.ArgumentTypeError(f"{x} not a valid angle, it should be between 0 and 359.")
    return x


def create_log(path, level):
    """Create the log as a text file and as a stream.

    :param path: the path of the log.
    :type path: str
    :param level: the level of the log.
    :type level: str
    :return: the logging:
    :rtype: logging
    """

    log_level_dict = {"DEBUG": logging.DEBUG,
                      "INFO": logging.INFO,
                      "WARNING": logging.WARNING,
                      "ERROR": logging.ERROR,
                      "CRITICAL": logging.CRITICAL}

    if level is None:
        log_level = log_level_dict["INFO"]
    else:
        log_level = log_level_dict[level]

    logging.basicConfig(format="%(asctime)s %(levelname)s:\t%(message)s",
                        datefmt="%Y/%m/%d %H:%M:%S",
                        level=log_level,
                        handlers=[logging.FileHandler(path), logging.StreamHandler()])
    return logging


def parse_frames(frames_selections, traj_files_paths):
    """
    Parse the frames' selection by trajectory file.

    :param frames_selections: the frames selection by trajectory file.
    :type frames_selections: str
    :param traj_files_paths: the trajectory files paths.
    :type traj_files_paths: list
    :return: the selected frames by trajectory file.
    :rtype: dict
    """
    frames_selection_data = {}
    # the pattern to get the beginning and the end of the frames selection
    pattern = re.compile("(.+):(\\d+|\\*)-(\\d+|\\*)")

    if frames_selections:
        start = 'the first frame'
        end = 'the last frame'
        traj_basenames = [os.path.basename(traj_files_path) for traj_files_path in traj_files_paths]
        for frames_sel in frames_selections.split(","):
            match = pattern.search(frames_sel)
            if match:
                current_traj = match.group(1)
                start = match.group(2)
                end = match.group(3)
                if current_traj not in traj_basenames:
                    raise argparse.ArgumentTypeError(f"The trajectory file {current_traj} in frame selection part "
                                                     f"'{frames_sel}' is not a file belonging to the inputs trajectory "
                                                     f"files: {','.join(traj_basenames)}")
                if start != "*":
                    frames_selection_data[current_traj] = {"begin": int(start)}
                if end != "*":
                    if current_traj not in frames_selection_data:
                        frames_selection_data[current_traj] = {"end": int(end)}
                    else:
                        frames_selection_data[current_traj]["end"] = int(end)
            else:
                raise argparse.ArgumentTypeError(f"The frame selection part '{frames_sel}' do not match the correct "
                                                 f"pattern {pattern.pattern}'")
        if comm.rank == 0:
            logging.info("Frames selection on input trajectory files:")
        for current_traj in frames_selection_data:
            logging.info(f"\t{current_traj}: frames selection from {start} to {end}.")
    return frames_selection_data


def resume_or_initialize_analysis(trajectory_files, topology_file, smp, distance_contacts, proportion_contacts,
                                  sim_time, resume_yaml, frames_sel):
    """
    Load the previous analysis data or create a new one if no previous analysis path was performed.

    :param trajectory_files: the current analysis trajectory files path.
    :type trajectory_files: list
    :param topology_file: the trajectories' topology file.
    :type topology_file: str
    :param smp: the sample name.
    :type smp: str
    :param distance_contacts: the threshold atoms distance in Angstroms for contacts.
    :type distance_contacts: float
    :param proportion_contacts: the minimal percentage of contacts for atoms contacts of different residues in the
    selected frames.
    :type proportion_contacts: float
    :param sim_time: the molecular dynamics simulation time in ns.
    :type sim_time: int
    :param resume_yaml: the path to the YAML file of previous analysis.
    :type resume_yaml: str
    :param frames_sel: the frames' selection for new trajectory files.
    :type frames_sel: dict
    :return: the initialized or resumed analysis data, the trajectory files' already analyzed.
    :rtype: dict, list
    """
    trajectory_files_to_skip = []
    if resume_yaml:
        with open(resume_yaml, "r") as file_handler:
            data = yaml.safe_load(file_handler.read())

        # check the processed analyzed trajectories files
        for t_file_path in trajectory_files:
            t_file = os.path.basename(t_file_path)
            if t_file in data["trajectory files processed"]:
                trajectory_files_to_skip.append(t_file)
                if comm.rank == 0:
                    logging.warning(f"\t{t_file} already processed in the previous analysis (check the YAML file, "
                                    f"section 'trajectory files processed'), this trajectory analysis is skipped.")

        discrepancies = []
        if data["sample"] != smp:
            discrepancies.append(f"discrepancy in --sample, current analysis is {smp}, previous analysis was "
                                 f"{data['sample']}")
        if data["topology file"] != os.path.basename(topology_file):
            discrepancies.append(f"discrepancy in --topology, current analysis is {topology_file}, previous analysis "
                                 f"was {data['topology file']}")
        if data["parameters"]["maximal atoms distance"] != distance_contacts:
            discrepancies.append(f"discrepancy in --distance-contacts, current analysis is {distance_contacts}, "
                                 f"previous analysis was {data['parameters']['maximal atoms distance']}")
        if data["parameters"]["proportion contacts"] != proportion_contacts:
            discrepancies.append(f"discrepancy in --proportion-contacts, current analysis is {proportion_contacts}, "
                                 f"previous analysis was {data['parameters']['proportion contacts']}")
        if discrepancies:
            discrepancies_txt = None
            for item in discrepancies:
                if discrepancies_txt:
                    discrepancies_txt = f"{discrepancies_txt}; {item}"
                else:
                    discrepancies_txt = item
                discrepancies_txt = f"{discrepancies_txt}. Check {resume_yaml}"
            raise KeyError(discrepancies_txt)

        # load the neighbors from the pickle file
        try:
            with open(data["pickle neighbors"], "rb") as file_handler:
                data["neighbors"] = pickle.load(file_handler)
        except FileNotFoundError as fnf_ex:
            logging.error(fnf_ex, exc_info=True)
            sys.exit(1)

        # add frames selection in new trajectory files
        if frames_sel:
            if "frames selections" in data["parameters"]:
                for traj_fn in frames_sel:
                    data["parameters"]["frames selections"][traj_fn] = frames_sel[traj_fn]
            else:
                data["parameters"]["frames selections"] = frames_sel
    else:
        data = {"sample": smp, "size Gb": 0, "frames": 0,
                "parameters": {"maximal atoms distance": distance_contacts, "proportion contacts": proportion_contacts},
                "topology file": os.path.basename(topology_file)}
    # set the simulation time
    data["parameters"]["time"] = f"{sim_time} ns"
    # add a neighbors section if necessary
    if "neighbors" not in data:
        data["neighbors"] = {}
    return data, trajectory_files_to_skip


def remove_processed_trajectories(all_traj, traj_to_skip, yaml_file):
    """
    Remove already processed trajectories.

    :param all_traj: the input trajectories.
    :type all_traj: list
    :param traj_to_skip: the trajectories already processed.
    :type traj_to_skip: list
    :param yaml_file: the path to the yaml file.
    :type yaml_file: str
    :return: the remaining trajectories.
    :rtype: list
    """
    traj_to_process = []
    for traj_path in all_traj:
        if os.path.basename(traj_path) not in traj_to_skip:
            traj_to_process.append(traj_path)
    if not traj_to_process:
        if comm.rank == 0:
            logging.warning(f"All trajectories ({','.join(all_traj)}) have already been processed, check the "
                            f"'trajectory files processed' section of: {yaml_file}")
    return traj_to_process


def load_trajectory(trajectory_file, topology_file, frames_sel):
    """
    Load the trajectory file.

    :param trajectory_file: the trajectory file path.
    :type trajectory_file: str
    :param topology_file: the topology file path.
    :type topology_file: str
    :param frames_sel: the frames' selection for new trajectory files.
    :type frames_sel: dict
    :return: the loaded trajectory.
    :rtype: pytraj.Trajectory
    """
    if comm.rank == 0:
        logging.info(f"\tLoading trajectory file, please be patient..")
    traj = None
    try:
        if trajectory_file in frames_sel:
            traj = pt.iterload(trajectory_file, top=topology_file,
                               frames_indices=range(frames_sel[trajectory_file]["begin"] - 1,
                                                    frames_sel[trajectory_file]["end"] - 1))
        else:
            traj = pt.iterload(trajectory_file, top=topology_file)
    except ValueError as ve_ex:
        logging.error(f"\tOne of the following files is missing: {trajectory_file} or {topology_file}")
        sys.exit(1)
    if comm.rank == 0:
        logging.info(f"\t\tMolecules:{traj.topology.n_mols:>20}")
        logging.info(f"\t\tResidues:{traj.topology.n_residues:>22}")
        logging.info(f"\t\tAtoms:{traj.topology.n_atoms:>27}")
        logging.info(f"\t\tTrajectory total frames:{traj.n_frames:>7}")
        logging.info(f"\t\tTrajectory memory size:{round(traj._estimated_GB, 6):>14} Gb")
    if os.path.basename(trajectory_file) in frames_sel:
        traj_bn = os.path.basename(trajectory_file)
        if frames_sel[traj_bn]["end"] > traj.n_frames:
            raise IndexError(f"Selected upper frame limit for {traj_bn} ({frames_sel[traj_bn]['end']}) from "
                             f"--frames argument is greater than the total frames number ({traj.n_frames}) of the MD "
                             f"trajectory.")
        frames_range = range(frames_sel[traj_bn]["begin"], frames_sel[traj_bn]["end"])
        if frames_sel[traj_bn]["begin"] == 1:
            frames_range[0] = 0
        traj = traj[frames_range]
        if comm.rank == 0:
            logging.info(f"\t\tSelected frames:{frames_sel[traj_bn]['begin']:>14} to {frames_sel[traj_bn]['end']}")
            logging.info(f"\t\tSelected frames memory size:{round(traj._estimated_GB, 6):>9} GB")
    else:
        txt = f"1 to {traj.n_frames}"
        if comm.rank == 0:
            logging.info(f"\t\tSelected frames:{txt:>20}")
    return traj


def check_trajectories_consistency(traj, path, data, frames_sel):
    """
    Check if the trajectory attributes match with the previous trajectories.

    :param traj: the current trajectory.
    :type traj: pytraj.Trajectory
    :param path: the current trajectory path.
    :type path: str
    :param data: the trajectories' data.
    :type data: dict
    :param frames_sel: the frames' selection on the trajectory files.
    :type frames_sel: dict
    :return: the updated trajectories' data.
    :rtype: dict
    """
    if "residues" not in data:
        data["residues"] = traj.topology.n_residues
        data["atoms"] = traj.topology.n_atoms
        data["molecules"] = traj.topology.n_mols
    else:
        if data["residues"] != traj.topology.n_residues:
            raise ValueError(f"the residues number ({traj.topology.n_residues}) is different from the residues number "
                             f"of the previous trajectories ({data['residues']}), check if "
                             f"{os.path.basename(path)} is from the same trajectory than the previous ones.")
        if data["atoms"] != traj.topology.n_atoms:
            raise ValueError(f"the atoms number ({traj.topology.n_atoms}) is different from the atoms number of "
                             f"the previous trajectories ({data['atoms']}), check if {os.path.basename(path)} is "
                             f"from the same trajectory than the previous ones.")
        if data["molecules"] != traj.topology.n_mols:
            raise ValueError(f"the molecules number ({traj.topology.n_mols}) is different from the molecules number of "
                             f"the previous trajectories ({data['molecules']}), check if {os.path.basename(path)} "
                             f"is from the same trajectory than the previous ones.")
    if frames_sel:
        if "frames selection" not in data:
            data["frames selection"] = {}
            for traj_name in frames_sel:
                data["frames selection"][traj_name] = frames_sel[traj_name]
    data["size Gb"] += traj._estimated_GB
    data["frames"] += traj.n_frames
    return data


def get_roi_atoms_indices(traj, roi):
    """
    get the region of interest atoms' indices.
    :param traj: the trajectory.
    :type traj: pytraj.Trajectory
    :param roi: the region of interest start and end residues index (1-index).
    :type roi: list
    :return: the region of interest atoms indices.
    :rtype: list
    """
    atoms_indices = []
    for index_atom in range(traj.topology.n_atoms):
        if traj.topology.atom(index_atom).resid in range(roi[0] - 1, roi[1]):
            atoms_indices.append(index_atom)
    return atoms_indices


def search_neighbors(roi_atoms_indices, traj, distance_thr, proportion_thr):
    """
    Get the region of interest atoms neighbors atoms.

    :param roi_atoms_indices: the region of interest atoms indices.
    :type roi_atoms_indices: list
    :param traj: the trajectory.
    :type traj: pytraj.Trajectory
    :param distance_thr: the maximal distance threshold in Angstroms to be considered as a neighbor.
    :type distance_thr: float
    :param proportion_thr: the frames proportion threshold to validate a contact.
    :type proportion_thr: float
    :return: the neighbors.
    :rtype: dict
    """
    data = {}
    for roi_atom_idx in roi_atoms_indices:
        data[roi_atom_idx] = {}
        mask = f"@{roi_atom_idx}<@{distance_thr}"
        neighbor_indices = pt.search_neighbors(traj, mask)
        for frame_indices in neighbor_indices:
            for atom_idx in list(frame_indices.values):
                if atom_idx in data[roi_atom_idx]:
                    data[roi_atom_idx][atom_idx] += 1
                else:
                    data[roi_atom_idx][atom_idx] = 1
        data[roi_atom_idx] = dict(sorted(data[roi_atom_idx].items()))
        threshold = traj.n_frames * (proportion_thr / 100.0)
        for atom_idx in list(data[roi_atom_idx].keys()):
            if data[roi_atom_idx][atom_idx] < threshold:
                del data[roi_atom_idx][atom_idx]
        if not data[roi_atom_idx]:
            del data[roi_atom_idx]
    return data


def get_neighbors_out_of_roi(data, atoms_neighbors, roi, padding, traj):
    """
    Get only the region of interest neighbors out of this region of interest.

    :param data: the neighbors' data.
    :type data: dict
    :param atoms_neighbors: the neighbors' atoms.
    :type atoms_neighbors: dict
    :param roi: the region of interest start and end residues index (1-index).
    :type roi: list
    :param padding: the number of residues to ignore at the region of interest borders.
    :type padding: int
    :param traj: the trajectory.
    :type traj: pytraj.Trajectory
    :return:
    """
    for roi_atom_idx in atoms_neighbors:
        idx_residue1 = traj.topology.atom(roi_atom_idx).resid
        residue1 = traj.topology.residue(idx_residue1).name
        atom1 = traj.topology.atom(roi_atom_idx).name
        for neighbor_atom_idx in atoms_neighbors[roi_atom_idx]:
            neighbor_atom = traj.topology.atom(neighbor_atom_idx)
            idx_residue2 = neighbor_atom.resid
            residue2 = traj.topology.residue(idx_residue2).name
            atom2 = traj.topology.atom(neighbor_atom_idx).name
            if idx_residue2 not in range(roi[0] - 1 - padding, roi[1] + padding):
                data["neighbors"].append(f"{residue1}{idx_residue1 + 1}_{atom1}-{residue2}{idx_residue2 + 1}_{atom2}")
                data["residue 1 position"].append(idx_residue1 + 1)
                data["residue 1"].append(residue1)
                data["atom 1"].append(atom1)
                data["residue 2 position"].append(idx_residue2 + 1)
                data["residue 2"].append(residue2)
                data["atom 2"].append(atom2)
                data["proportion frames (%)"].append(
                    round(atoms_neighbors[roi_atom_idx][neighbor_atom_idx] / traj.n_frames * 100, 1))

    pd.set_option('display.max_rows', None)
    return pd.DataFrame.from_dict(data)


def record_analysis(data, out_dir, current_trajectory_file, smp):
    """
    Record the analysis to a YAML file and pickle the neighbors' data in a binary file.

    :param data: the trajectory analysis data.
    :type data: dict
    :param out_dir: the path to the output directory.
    :type out_dir: str
    :param current_trajectory_file: the current trajectory file.
    :type current_trajectory_file: str
    :param smp: the sample name.
    :type smp: str
    :return: the trajectory analysis data.
    :rtype: dict
    """
    if "trajectory files processed" not in data:
        data["trajectory files processed"] = []
    data["trajectory files processed"].append(os.path.basename(current_trajectory_file))

    # extract the neighbors from the data and pickle them
    out_pickle = os.path.join(out_dir, f"{smp.replace(' ', '_')}_analysis.pkl")
    neighbors_analysis_data = data.pop("neighbors")
    with open(out_pickle, "wb") as file_handler:
        pickle.dump(neighbors_analysis_data, file_handler)
    if comm.rank == 0:
        logging.info(f"\t\tNeighbors analysis saved: {os.path.abspath(out_pickle)}")

    # save the analysis parameters without the neighbors to a YAML file
    data["pickle neighbors"] = os.path.abspath(out_pickle)
    out_yaml = os.path.join(out_dir, f"{smp.replace(' ', '_')}_analysis.yaml")
    with open(out_yaml, "w") as file_handler:
        yaml.dump(data, file_handler)
    if comm.rank == 0:
        logging.info(f"\t\tAnalysis parameters saved: {os.path.abspath(out_yaml)}")

    # add the extracted neighbors to the data again
    data["neighbors"] = neighbors_analysis_data
    return data


if __name__ == "__main__":
    descr = f"""
    {os.path.basename(__file__)} v. {__version__}

    Created by {__author__}.
    Contact: {__email__}
    {__copyright__}

    Distributed on an "AS IS" basis without warranties or conditions of any kind, either express or implied.

    From molecular dynamics trajectories files (*.nc), the script performs a trajectory analysis to search contacts 
    between atoms.
    """
    parser = argparse.ArgumentParser(description=descr, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--out", required=True, type=str, help="the path to the output directory.")
    parser.add_argument("-s", "--sample", required=True, type=str,
                        help="the sample ID used for the files name.")
    parser.add_argument("-t", "--topology", required=True, type=str,
                        help="the path to the molecular dynamics topology file.")
    parser.add_argument("-n", "--nanoseconds", required=True, type=int,
                        help="the molecular dynamics simulation time in nano seconds.")
    parser.add_argument("-r", "--roi", required=True, type=int, nargs=2,
                        help="the region of interest start and end residues positions (i-index). The atoms of this "
                             "region will be excluded from the contacts search.")
    parser.add_argument("-x", "--padding", required=False, type=int, default=5,
                        help="the minimal number of residues to exclude from the start and end residues positions to "
                             "search. If the start position is 20 and the end position is 37 and the padding is 5, the "
                             "neighbors search will be between the 20th and the 37th residues and the 1st to 15th "
                             "residues and the 20th and the 37th residues and the 42th and last residue.")
    parser.add_argument("-f", "--frames", required=False, type=str,
                        help="the frames selection by trajectory file. The arguments should be <TRAJ_FILE>:100-1000. "
                             "If the <TRAJ_FILE> contains 2000 frames, only the frames from 100-1000 will be selected."
                             "Multiple frames selections can be performed with comma separators, i.e: "
                             "'<TRAJ_FILE_1>:100-1000,<TRAJ_FILE_2>:1-500'. If a '*' is used as '<TRAJ_FILE_1>:100-*', "
                             "the used frames will be the 100th until the end frame.")
    parser.add_argument("-d", "--distance-contacts", required=False, type=restricted_positive, default=3.0,
                        help="The minimal contact distance between two atoms. Default is 3.0 Angstroms.")
    parser.add_argument("-p", "--proportion-contacts", required=False, type=restricted_float, default=20.0,
                        help="the minimal percentage of frames which make contact between 2 atoms of different "
                             "residues in the selected frame of the molecular dynamics simulation, default is 20%%.")
    parser.add_argument("--resume", required=False, type=str,
                        help="the YAML file path of the previous trajectory analysis. The analysis of the new "
                             "trajectory files of the same system will resume on the previous trajectory analysis. The "
                             "new trajectory file analysis will be added to the previous data.")
    parser.add_argument("-l", "--log", required=False, type=str,
                        help="the path for the log file. If this option is skipped, the log file is created in the "
                             "output directory.")
    parser.add_argument("--log-level", required=False, type=str,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="set the log level. If the option is skipped, log level is INFO.")
    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument("inputs", nargs="+", type=str,
                        help="the paths to the molecular dynamics trajectory files.")
    args = parser.parse_args()

    comm = MPI.COMM_WORLD
    if comm.rank == 0:
        # create output directory if necessary
        os.makedirs(args.out, exist_ok=True)

        # create the logger
        if args.log:
            log_path = args.log
        else:
            log_path = os.path.join(args.out, f"{os.path.splitext(os.path.basename(__file__))[0]}.log")
        create_log(log_path, args.log_level)

        logging.info(f"version: {__version__}")
        logging.info(f"CMD: {' '.join(sys.argv)}")
        logging.info(f"Atoms neighboring search for the residues region of interest: {args.roi[0]:>3}-{args.roi[1]}")
        logging.info(f"Atoms maximal contacts distance threshold: {args.distance_contacts:>24} \u212B")
        logging.info(f"Applied mask for each \"X\" region of interest's atom number:\t@X<@{args.distance_contacts}")
        logging.info(f"Minimal frames proportion with atoms contacts: {args.proportion_contacts:>21.1f}%")
        logging.info(f"Molecular Dynamics duration: {args.nanoseconds:>36} ns")

    frames_selection = None
    try:
        frames_selection = parse_frames(args.frames, args.inputs)
    except argparse.ArgumentTypeError as ex:
        logging.error(ex, exc_info=True)
        sys.exit(1)

    data_traj = None
    traj_files = None
    try:
        data_traj, skipped_traj = resume_or_initialize_analysis(args.inputs, args.topology, args.sample,
                                                                args.distance_contacts, args.proportion_contacts,
                                                                args.nanoseconds, args.resume, frames_selection)
        traj_files = remove_processed_trajectories(args.inputs, skipped_traj, args.resume)
    except KeyError as exc:
        logging.error(exc, exc_info=True)
        sys.exit(1)

    trajectory = None
    neighbors = {"neighbors": [], "residue 1 position": [], "residue 1": [], "atom 1": [], "residue 2 position": [],
                 "residue 2": [], "atom 2": [], "proportion frames (%)": []}
    for traj_file in traj_files:
        # load the trajectory
        if comm.rank == 0:
            logging.info(f"Processing trajectory file: {traj_file}")
        try:
            trajectory = load_trajectory(traj_file, args.topology, frames_selection)
            data_traj = check_trajectories_consistency(trajectory, traj_file, data_traj, frames_selection)
        except RuntimeError as exc:
            logging.error(f"Check if the topology ({args.topology}) and/or the trajectory ({', '.join(args.inputs)}) "
                          f"files exists", exc_info=True)
            sys.exit(1)
        except ValueError as exc:
            logging.error(exc, exc_info=True)
            sys.exit(1)
        except IndexError as exc:
            logging.error(exc, exc_info=True)
            sys.exit(1)

        # find neighbors
        logging.info(f"\tSearching atoms neighbors of the residues region of interest ({args.roi[0]}-{args.roi[1]}) "
                     f"atoms with {args.padding} residue{'s' if args.padding > 1 else ''} padding:")

        atoms_indices_in_roi = get_roi_atoms_indices(trajectory, args.roi)
        atoms_neighbors_with_roi_atoms = search_neighbors(atoms_indices_in_roi, trajectory, args.distance_contacts,
                                                          args.proportion_contacts)
        neighbors = get_neighbors_out_of_roi(neighbors, atoms_neighbors_with_roi_atoms, args.roi, args.padding, trajectory)
        logging.info(f"\t\tAtoms neighbor{'s' if len(neighbors) > 1 else ''} found:{len(neighbors):>15}")

        # pickle the analysis
        record_analysis(data_traj, args.out, traj_file, args.sample)

    if comm.rank == 0:
        logging.info(f"{len(data_traj['trajectory files processed'])} processed trajectory files: "
                     f"{', '.join(data_traj['trajectory files processed'])}")
        logging.info(f"Whole trajectories memory size: {round(data_traj['size Gb'], 6):>17} Gb")
        logging.info(f"Whole trajectories frames: {data_traj['frames']:>16}")
        logging.info(f"Whole trajectories neighbors found: {len(neighbors)}")
        out_path = os.path.join(args.out, f"neighbors_{args.sample.replace(' ', '_')}.csv")
        neighbors.to_csv(out_path, index=False)
        logging.info(f"Neighbors result file: {os.path.abspath(out_path)}")

