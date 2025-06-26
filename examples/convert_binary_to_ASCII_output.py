#!/usr/bin/env python3
"""
Binary to ASCII Converter for JETSCAPE Final State Files

This script reads a binary file containing event data in a specific format and writes
an ASCII representation of the data to a specified output file. It supports both
parton and hadron final states and is suitable for large files through buffered streaming.

Usage:
    python convert_binary_to_ascii.py input_file output_file partons/hadrons

Author: Hendrik Roch
"""

import argparse
import struct


def convert_binary_to_ASCII(
    binary_file_path, ascii_file_path, parton_or_hadron
):
    """
    Convert a binary final state file to ASCII format.

    Parameters
    ----------
    binary_file_path : str
        Path to the binary file.
    ascii_file_path : str
        Output path for the ASCII file.
    parton_or_hadron : str
        Either 'partons' or 'hadrons' indicating the final state type.

    Raises
    ------
    ValueError
        If the file is too short or malformed.
    """
    with open(binary_file_path, "rb") as bin_f, open(
        ascii_file_path, "w", buffering=1024 * 1024
    ) as out_f:
        # --- HEADER ---
        raw = bin_f.read(20)
        if len(raw) < 20:
            raise ValueError("File too short to contain valid header.")

        file_version, *status_codes = struct.unpack("5i", raw)
        out_f.write(
            f"JETSCAPE_FINAL_STATE\tv{file_version}\t|\tN\tpid\tstatus\tE\tPx\tPy\tPz\n"
        )
        # Skipped status codes are not written to the output file, write only message to terminal
        # Only print if status codes are different from 666
        for code in status_codes:
            if code != 666:
                print(f"Status code {code} was filtered out in original file.")

        # --- Events ---
        bin_f.seek(0, 2)
        file_size = bin_f.tell()
        FOOTER_SIZE = 16  # two doubles
        bin_f.seek(20)  # After header
        offset = 20

        while offset < file_size - FOOTER_SIZE:
            bin_f.seek(offset)
            header = bin_f.read(36)
            if len(header) < 36:
                raise ValueError(
                    f"Unexpected EOF reading event header at offset {offset}"
                )

            (event_number,) = struct.unpack_from("i", header, 0)
            (event_weight,) = struct.unpack_from("d", header, 4)
            event_plane_angle, vx, vy, vz, centrality, pt_hat = (
                struct.unpack_from("6f", header, 12)
            )

            offset += 36
            bin_f.seek(offset)
            raw = bin_f.read(4)
            if len(raw) < 4:
                raise ValueError(
                    f"Unexpected EOF reading particle count at offset {offset}"
                )
            (num_particles,) = struct.unpack("i", raw)
            offset += 4

            particle_data_size = 24 * num_particles
            bin_f.seek(offset)
            particle_data = bin_f.read(particle_data_size)
            if len(particle_data) < particle_data_size:
                raise ValueError(
                    f"Unexpected EOF reading {num_particles} particles at offset {offset}"
                )
            offset += particle_data_size

            out_f.write(
                f"#\tEvent\t{event_number}\t"
                f"weight\t{event_weight:.15f}\t"
                f"EPangle\t{event_plane_angle:.6f}\t"
                f"N_{parton_or_hadron}\t{num_particles}\t"
                f"vertex_x\t{vx:.6f}\tvertex_y\t{vy:.6f}\tvertex_z\t{vz:.6f}\t"
                f"centrality\t{centrality:.6f}\tpt_hat\t{pt_hat:.6f}\n"
            )

            for i in range(num_particles):
                base = i * 24
                pid, status = struct.unpack_from("2i", particle_data, base)
                e, px, py, pz = struct.unpack_from(
                    "4f", particle_data, base + 8
                )
                out_f.write(
                    f"{i} {pid} {status} "
                    f"{e:.6f} {px:.6f} {py:.6f} {pz:.6f}\n"
                )

        # --- Footer ---
        bin_f.seek(offset)
        footer = bin_f.read(16)
        if len(footer) != 16:
            raise ValueError("Unexpected EOF reading footer.")
        sigmaGen, sigmaErr = struct.unpack("2d", footer)
        out_f.write(f"#\tsigmaGen\t{sigmaGen:.6e}\tsigmaErr\t{sigmaErr:.6e}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert a JETSCAPE binary final state file to ASCII format."
    )
    parser.add_argument(
        "binary_file_path", type=str, help="Path to the binary file to convert."
    )
    parser.add_argument(
        "ascii_file_path",
        type=str,
        help="Path to save the converted ASCII file.",
    )
    parser.add_argument(
        "parton_or_hadron",
        type=str,
        choices=["partons", "hadrons"],
        help="Specify 'partons' or 'hadrons' for the conversion.",
    )
    args = parser.parse_args()
    # Call the function to convert the binary file to ASCII
    convert_binary_to_ASCII(
        args.binary_file_path, args.ascii_file_path, args.parton_or_hadron
    )
    # Print a message indicating the conversion is complete
    print(
        f"Conversion complete. The ASCII file is saved at {args.ascii_file_path}."
    )
