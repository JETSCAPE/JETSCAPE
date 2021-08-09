#!/usr/bin/env python3

""" Parse JETSCAPE ascii input files in chunks.

.. codeauthor:: Raymond Ehlers
"""

import itertools
import logging
import os
import typing
from pathlib import Path
from typing import Any, Callable, Iterator, List, Optional, Union

import awkward as ak
import attr
import numpy as np


logger = logging.getLogger(__name__)


class ReachedEndOfFileException(Exception):
    """Indicates that we've somehow hit the end of the file.

    We have a separate exception so we can pass additional information
    about the context if desired.
    """
    ...


class ReachedXSecAtEndOfFileException(ReachedEndOfFileException):
    """Indicates that we've hit the cross section in the last line of the file."""
    ...


@attr.s(frozen=True)
class CrossSection:
    value: float = attr.ib()
    error: float = attr.ib()


@attr.s(frozen=True)
class HeaderInfo:
    event_number: int = attr.ib()
    event_plane_angle: float = attr.ib()
    n_particles: int = attr.ib()
    event_weight: float = attr.ib(default=-1)


def _retrieve_last_line_of_file(f: typing.TextIO, read_chunk_size: int = 100) -> str:
    """ Retrieve the last line of the file.

    From: https://stackoverflow.com/a/7167316/12907985

    Args:
        f: File-like object.
        read_chunk_size: Size of step in bytes to read backwards into the file. Default: 100.

    Returns:
        Last line of file, assuming it's found.
    """
    last_line = ""
    while True:
        # We grab chunks from the end of the file towards the beginning until we
        # get a new line
        # However, we need to be more careful because seeking directly from the end
        # of a text file is apparently undefined behavior. To address this, we
        # follow the approach from here: https://stackoverflow.com/a/51131242/12907985
        # First, move to the end of the file.
        f.seek(0, os.SEEK_END)
        # Then, just back from the current position (based on SEEK_SET and tell()
        # of the current position).
        f.seek(f.tell() - len(last_line) - read_chunk_size, os.SEEK_SET)
        # NOTE: This chunk isn't necessarily going back read_chunk_size characters, but
        #       rather just some number of bytes. Depending on the encoding, it may be.
        #       In any case, we don't care - we're just looking for the last line.
        chunk = f.read(read_chunk_size)

        if not chunk:
            # The whole file is one big line
            return last_line

        if not last_line and chunk.endswith('\n'):
            # Ignore the trailing newline at the end of the file (but include it
            # in the output).
            last_line = '\n'
            chunk = chunk[:-1]

        nl_pos = chunk.rfind('\n')
        # What's being searched for will have to be modified if you are searching
        # files with non-unix line endings.

        last_line = chunk[nl_pos + 1:] + last_line

        if nl_pos == -1:
            # The whole chunk is part of the last line.
            continue

        return last_line


def _parse_cross_section(line: str) -> CrossSection:
    # Parse the string.
    values = line.split("\t")
    if len(values) == 5 and values[1] == "sigmaGen":
        ###################################
        # Cross section info specification:
        ###################################
        # The cross section info is formatted as follows, with each entry separated by a `\t` character:
        # # sigmaGen 182.423 sigmaErr 11.234
        # 0 1        2       3        4
        #
        # I _think_ the units are mb^-1
        info = CrossSection(
            value=float(values[2]),     # Cross section
            error=float(values[4]),     # Cross section error
        )
    else:
        raise ValueError(f"Parsing of comment line failed: {values}")

    return info


def _extract_x_sec_and_error(f: typing.TextIO, read_chunk_size: int = 100) -> Optional[CrossSection]:
    """ Extract cross section and error from the end of the file.

    Args:
        f: File-like object.
        read_chunk_size: Size of step in bytes to read backwards into the file. Default: 100.

    Returns:
        Cross section and error, if found.
    """
    # Retrieve the last line of the file.
    last_line = _retrieve_last_line_of_file(f=f, read_chunk_size=read_chunk_size)
    # Move the file back to the start to reset it for later use.
    f.seek(0)

    # logger.debug(f"last line: {last_line}")
    if last_line.startswith("#\tsigmaGen"):
        # logger.debug("Parsing xsec")
        return _parse_cross_section(line = last_line)

    return None


def _parse_header_line_format_unspecified(line: str) -> HeaderInfo:
    """Parse line that is expected to be a header.

    The most common case is that it's a header, in which case we parse the line. If it's not a header,
    we also check if it's a cross section, which we note as having found at the end of the file via
    an exception.

    Args:
        line: Line to be parsed.
    Returns:
        The HeaderInfo information that was extracted.
    Raises:
        ReachedXSecAtEndOfFileException: If we find the cross section.
    """
    # Parse the string.
    values = line.split("\t")
    # Compare by length first so we can short circuit immediately if it doesn't match, which should
    # save some string comparisons.
    info: Union[HeaderInfo, CrossSection]
    if (len(values) == 19 and values[1] == "Event") or (len(values) == 17 and values[1] == "Event"):
        ##########################
        # Header v2 specification:
        ##########################
        # As of 20 April 2021, the formatting of the header has been improved.
        # This function was developed to parse it.
        # The header is defined as follows, with each entry separated by a `\t` character:
        # `# Event 1 weight 0.129547 EPangle 0.0116446 N_hadrons 236 | N  pid status E Px Py Pz Eta Phi`
        #  0 1     2 3      4        5       6         7         8   9 10 11  12    13 14 15 16 17  18
        #
        # NOTE: Everything after the "|" is just documentation for the particle entries stored below.
        #
        info = HeaderInfo(
            event_number=int(values[2]),            # Event number
            event_plane_angle=float(values[6]),     # EP angle
            n_particles=int(values[8]),             # Number of particles
            event_weight=float(values[4]),          # Event weight
        )
    elif len(values) == 9 and "Event" in values[2]:
        ##########################
        # Header v1 specification:
        ##########################
        # The header v1 specification is as follows, with ">" followed by same spaces indicating a "\t" character:
        # #>  0.0116446>  Event1ID>   236>pstat-EPx>  Py> Pz> Eta>Phi
        # 0   1           2           3   4           5   6   7   8
        #
        # NOTE: There are some difficulties in parsing this header due to inconsistent spacing.
        #
        # The values that we want to extract are:
        # index: Meaning
        #   1: Event plane angle. float, potentially in scientific notation.
        #   2: Event number, of the from "EventNNNNID". Can be parsed as val[5:-2] to generically extract `NNNN`. int.
        #   3: Number of particles. int.
        info = HeaderInfo(
            event_number=int(values[2][5:-2]),      # Event number
            event_plane_angle=float(values[1]),     # EP angle
            n_particles=int(values[3]),             # Number of particles
        )
    elif len(values) == 5 and values[1] == "sigmaGen":
        # If we've hit the cross section, and we're not doing the initial extraction of the cross
        # section, this means that we're at the last line of the file, and should notify as such.
        # NOTE: By raising with the cross section, we make it possible to retrieve it, even though
        #       we've raised an exception here.
        raise ReachedXSecAtEndOfFileException(_parse_cross_section(line))
    else:
        raise ValueError(f"Parsing of comment line failed: {values}")

    return info


def _parse_header_line_format_v2(line: str) -> HeaderInfo:
    """Parse line that is expected to be a header according to the v2 file format.

    The most common case is that it's a header, in which case we parse the line. If it's not a header,
    we also check if it's a cross section, which we note as having found at the end of the file via
    an exception.

    Args:
        line: Line to be parsed.
    Returns:
        The HeaderInfo information that was extracted.
    Raises:
        ReachedXSecAtEndOfFileException: If we find the cross section.
    """
    # Parse the string.
    values = line.split("\t")
    # Compare by length first so we can short circuit immediately if it doesn't match, which should
    # save some string comparisons.
    info: Union[HeaderInfo, CrossSection]
    if len(values) == 9 and values[1] == "Event":
        ##########################
        # Header v2 specification:
        ##########################
        # As of 22 June 2021, the formatting of the header is as follows:
        # This function was developed to parse it.
        # The header is defined as follows, with each entry separated by a `\t` character:
        # `# Event 1 weight 0.129547 EPangle 0.0116446 N_hadrons 236`
        #  0 1     2 3      4        5       6         7         8
        #
        info = HeaderInfo(
            event_number=int(values[2]),            # Event number
            event_plane_angle=float(values[6]),     # EP angle
            n_particles=int(values[8]),             # Number of particles
            event_weight=float(values[4]),          # Event weight
        )
    elif len(values) == 5 and values[1] == "sigmaGen":
        # If we've hit the cross section, and we're not doing the initial extraction of the cross
        # section, this means that we're at the last line of the file, and should notify as such.
        # NOTE: By raising with the cross section, we make it possible to retrieve it, even though
        #       we've raised an exception here.
        raise ReachedXSecAtEndOfFileException(_parse_cross_section(line))
    else:
        raise ValueError(f"Parsing of comment line failed: {values}")

    return info


# Register header parsing functions
_file_format_version_to_header_parser = {
    2: _parse_header_line_format_v2,
    -1: _parse_header_line_format_unspecified,
}


def _parse_event(f: Iterator[str], parse_header_line: Callable[[str], HeaderInfo]) -> Iterator[Union[HeaderInfo, str]]:
    """Parse a single event in a FinalState* file.

    Raises:
        ReachedXSecAtEndOfFileException: We've found the line with the xsec and error at the
            end of the file. Effectively, we've exhausted the iterator.
        ReachedEndOfFileException: We've hit the end of file without finding the xsec and
            error. This may be totally fine, depending on the version of the FinalState*
            output.
    """
    # Our first line should be the header, which will be denoted by a "#".
    # Let the calling know if there are no events left due to exhausting the iterator.
    try:
        header = parse_header_line(next(f))
        # logger.info(f"header: {header}")
        yield header
    except StopIteration:
        # logger.debug("Hit end of file exception!")
        raise ReachedEndOfFileException()

    # From the header, we know how many particles we have in the event, so we can
    # immediately yield the next n_particles lines. And this will leave the generator
    # at the next event, or (if we've exhausted the iterator) at the end of file (either
    # at the xsec and error, or truly exhausted).
    yield from itertools.islice(f, header.n_particles)


class ChunkNotReadyException(Exception):
    """Indicate that the chunk hasn't been parsed yet, and therefore is not ready."""
    ...


@attr.s
class ChunkGenerator:
    """ Generator a chunk of the file.

    Args:
        g: Iterator over the input file.
        events_per_chunk: Number of events for the chunk.
        cross_section: Cross section information.
        file_format_version: File format version. Default: -1, which corresponds to before the format
            was defined, and it will try it's best to guess the format.
    """
    g: Iterator[str] = attr.ib()
    _events_per_chunk: int = attr.ib()
    cross_section: Optional[CrossSection] = attr.ib(default=None)
    _file_format_version: int = attr.ib(default=-1)
    _headers: List[HeaderInfo] = attr.ib(factory=list)
    _reached_end_of_file: bool = attr.ib(default=False)

    def _is_chunk_ready(self) -> bool:
        """True if the chunk is ready"""
        return not (len(self._headers) == 0 and not self._reached_end_of_file)

    def _require_chunk_ready(self) -> None:
        """Require that the chunk is ready (ie. been parsed).

        Raises:
            ChunkNotReadyException: Raised if the chunk isn't ready.
        """
        if not self._is_chunk_ready():
            raise ChunkNotReadyException()

    @property
    def events_per_chunk(self) -> int:
        return self._events_per_chunk

    @property
    def reached_end_of_file(self) -> bool:
        self._require_chunk_ready()
        return self._reached_end_of_file

    @property
    def events_contained_in_chunk(self) -> int:
        self._require_chunk_ready()
        return len(self._headers)

    @property
    def headers(self) -> List[HeaderInfo]:
        self._require_chunk_ready()
        return self._headers

    def n_particles_per_event(self) -> np.ndarray:
        self._require_chunk_ready()
        return np.array([
            header.n_particles for header in self._headers
        ])

    def event_split_index(self) -> np.ndarray:
        self._require_chunk_ready()
        # NOTE: We skip the last header due to the way that np.split works.
        #       It will go from the last index to the end of the array.
        return np.cumsum([
            header.n_particles for header in self._headers
        ])[:-1]

    @property
    def incomplete_chunk(self) -> bool:
        self._require_chunk_ready()
        return len(self._headers) != self._events_per_chunk

    def __iter__(self) -> Iterator[str]:
        _parse_header_line = _file_format_version_to_header_parser[self._file_format_version]
        for _ in range(self._events_per_chunk):
            # logger.debug(f"i: {i}")
            try:
                event_iter = _parse_event(
                    self.g,
                    parse_header_line=_parse_header_line,
                )
                # NOTE: Typing gets ignored here because I don't know how to encode the additional
                #       information that the first yielded line will be the header, and then the rest
                #       will be strings. So we just ignore typing here. In principle, we could split
                #       up the _parse_event function, but I find condensing it into a single function
                #       to be more straightforward from a user perspective.
                # First, get the header. We know this first line must be a header
                self._headers.append(
                    next(event_iter)  # type: ignore
                )
                # Then we yield the rest of the particles in the event
                yield from event_iter  # type: ignore
            except (ReachedEndOfFileException, ReachedXSecAtEndOfFileException):
                # If we're reached the end of file, we should note that inside the chunk
                # because it may not have reached the full set of events per chunk.
                self._reached_end_of_file = True
                # Since we've reached this point, we need to stop.
                break
            # NOTE: If we somehow reached StopIteration, it's also fine - just
            #       allow it to propagate through and end the for loop.


def read_events_in_chunks(filename: Path, events_per_chunk: int = int(1e5)) -> Iterator[ChunkGenerator]:
    """ Read events in chunks from stored JETSCAPE FinalState* ASCII files.

    This provides access to the lines of the file itself, but it is up to the user to parse each line.
    Consequently, many useful features are implemented on top of it. Users are encouraged to use those
    more full featured functions, such as `read(...)`.

    Args:
        filename: Path to the file.
        events_per_chunk: Number of events to store in each chunk. Default: 1e5.
    Returns:
        Chunks iterator. When this iterator is consumed, it will generate lines from the file until it
            hits the number of events mark. The header information is contained inside the object.
    """
    # Validation
    filename = Path(filename)

    with open(filename, "r") as f:
        # First step, extract the final cross section and header.
        cross_section = _extract_x_sec_and_error(f)

        # Define an iterator so we can increment it in different locations in the code.
        # Fine to use if it the entire file fits in memory.
        #read_lines = iter(f.readlines())
        # Use this if the file doesn't fit in memory (fairly likely for these type of files)
        read_lines = iter(f)

        # Check for file format version indcating how we should parse it.
        file_format_version = -1
        first_line = next(read_lines)
        first_line_split = first_line.split("\t")
        if len(first_line_split) > 3 and first_line_split[1] == "JETSCAPE_FINAL_STATE":
            # 1: is to remove the "v" in the version
            file_format_version = int(first_line_split[2][1:])
        else:
            # We need to move back to the beginning of the file, so we just burned through
            # a meaningful line (which almost certianly contains an event header).
            # NOTE: My initial version of two separate iterators doesn't work because it appears
            #       that you cannot do so for a file (which I suppose I can make sense of because
            #       it points to a position in a file, but still unexpected).
            f.seek(0)

        logger.info(f"Found file format version: {file_format_version}")

        # Now, need to setup chunks.
        # NOTE: The headers and additional info are passed through the ChunkGenerator.
        while True:
            # We keep an explicit reference to the chunk so we can set the end of file state
            # if we reached the end of the file.
            chunk = ChunkGenerator(
                g=read_lines,
                events_per_chunk=events_per_chunk,
                cross_section=cross_section,
                file_format_version=file_format_version,
            )
            yield chunk
            if chunk.reached_end_of_file:
                break


class FileLikeGenerator:
    """ Wrapper class to make a generator look like a file.

    Pandas requires passing a filename or a file-like object, but we handle the generator externally
    so we can find each chunk boundary, parse the headers, etc. Consequently, we need to make this
    generator appear as if it's a file.

    Based on https://stackoverflow.com/a/18916457/12907985

    Args:
        g: Generator to be wrapped.
    """
    def __init__(self, g: Iterator[str]):
        self.g = g

    def read(self, n: int = 0) -> Any:
        """ Read method is required by pandas. """
        try:
            return next(self.g)
        except StopIteration:
            return ''

    def __iter__(self) -> Iterator[str]:
        """ Iteration is required by pandas. """
        return self.g


def _parse_with_pandas(chunk_generator: Iterator[str]) -> np.ndarray:
    """ Parse the lines with `pandas.read_csv`

    `read_csv` uses a compiled c parser. As of 6 October 2020, it is tested to be the fastest option.

    Args:
        chunk_generator: Generator of chunks of the input file for parsing.
    Returns:
        Array of the particles.
    """
    # Delayed import so we only take the import time if necessary.
    import pandas as pd

    return pd.read_csv(
        FileLikeGenerator(chunk_generator),
        # NOTE: If the field is missing (such as eta and phi), they will exist, but they willl be filled with NaN
        #       We actively take advantage of this so we don't have to change the parsing for header v1 (which
        #       includes eta and phi) vs header v2 (which does not)
        names=["particle_index", "particle_ID", "status", "E", "px", "py", "pz", "eta", "phi"],
        header=None,
        comment="#",
        sep=r"\s+",
        # Converting to numpy makes the dtype conversion moot.
        # dtype={
        #     "particle_index": np.int32, "particle_ID": np.int32, "status": np.int8,
        #     "E": np.float32, "px": np.float32, "py": np.float32, "pz": np.float32,
        #     "eta": np.float32, "phi": np.float32
        # },
        # We can reduce the number of columns when reading.
        # However, it makes little difference, makes it less general, and we can always drop the columns later.
        # So we disable it for now.
        # usecols=["particle_ID", "status", "E", "px", "py", "eta", "phi"],
        #
        # NOTE: It's important that we convert to numpy before splitting. Otherwise, it will return columns names,
        #       which will break the header indexing and therefore the conversion to awkward.
    ).to_numpy()


def _parse_with_python(chunk_generator: Iterator[str]) -> np.ndarray:
    """ Parse the lines with python.

    We have this as an option because np.loadtxt is surprisingly slow.

    Args:
        chunk_generator: Generator of chunks of the input file for parsing.
    Returns:
        Array of the particles.
    """
    particles = []
    for p in chunk_generator:
        if not p.startswith("#"):
            particles.append(np.array(p.rstrip("\n").split(), dtype=np.float64))
    return np.stack(particles)


def _parse_with_numpy(chunk_generator: Iterator[str]) -> np.ndarray:
    """ Parse the lines with numpy.

    Unfortunately, this option is surprisingly, presumably because it has so many options.
    Pure python appears to be about 2x faster. So we keep this as an option for the future,
    but it is not used by default.

    Args:
        chunk_generator: Generator of chunks of the input file for parsing.
    Returns:
        Array of the particles.
    """
    return np.loadtxt(chunk_generator)


def read(filename: Union[Path, str], events_per_chunk: int, parser: str = "pandas") -> Iterator[ak.Array]:
    """ Read a JETSCAPE FinalState{Hadrons,Partons} ASCII output file in chunks.

    This is the primary user function. We read in chunks to keep the memory usage manageable.

    Note:
        We store the data in the smallest possible types that can still encompass their range.

    Args:
        filename: Filename of the ASCII file.
        events_per_chunk: Number of events to provide in each chunk.
        parser: Name of the parser to use. Default: `pandas`, which uses `pandas.read_csv`. It uses
            compiled c, and seems to be the fastest available option. Other options: ["python", "numpy"].
    Returns:
        Generator of an array of events_per_chunk events.
    """
    # Validation
    filename = Path(filename)

    # Setup
    parsing_function_map = {
        "pandas": _parse_with_pandas,
        "python": _parse_with_python,
        "numpy": _parse_with_numpy,
    }
    parsing_function = parsing_function_map[parser]

    # Read the file, creating chunks of events.
    for i, chunk_generator in enumerate(read_events_in_chunks(filename=filename, events_per_chunk=events_per_chunk)):
        # Give a notification just in case the parsing is slow...
        logger.debug(f"New chunk {i}")

        # First, parse the lines. We need to make this call before attempt to convert into events because the necessary
        # info (namely, n particles per event) is only available and valid after we've parse the lines.
        res = parsing_function(iter(chunk_generator))

        # Before we do anything else, if our events_per_chunk is a even divisor of the total number of events
        # and we've hit the end of the file, we can return an empty generator after trying to parse the chunk.
        # In that case, we're done - just break.
        # NOTE: In the case where we've reached the end of file, but we've parsed some events, we want to continue
        #       on so we don't lose any events.
        if chunk_generator.reached_end_of_file and len(chunk_generator.headers) == 0:
            break

        # Now, convert into the awkward array structure.
        array_with_events = ak.unflatten(
            ak.Array(res), chunk_generator.n_particles_per_event()
        )

        # Cross checks.
        # Length check that we have as many events as expected based on the number of headers.
        # logger.debug(f"ak.num: {ak.num(array_with_events, axis = 0)}, len headers: {len(chunk_generator.headers)}")
        assert (ak.num(array_with_events, axis = 0) == len(chunk_generator.headers))
        # Check that n particles agree
        n_particles_from_header = np.array([header.n_particles for header in chunk_generator.headers])
        # logger.info(f"n_particles from headers: {n_particles_from_header}")
        # logger.info(f"n_particles from array: {ak.num(array_with_events, axis = 1)}")
        assert (np.asarray(ak.num(array_with_events, axis = 1)) == n_particles_from_header).all()
        # State of the chunk
        # logger.debug(f"Reached end of file: {chunk_generator.reached_end_of_file}")
        # logger.debug(f"Incomplete chunk: {chunk_generator.incomplete_chunk}")
        # Let the use know so they're not surprised.
        if chunk_generator.incomplete_chunk:
            logger.warning(f"Requested {chunk_generator.events_per_chunk} events, but only {chunk_generator.events_contained_in_chunk} are available because we hit the end of the file.")

        # Header info
        header_level_info = {
            "event_plane_angle": np.array([header.event_plane_angle for header in chunk_generator.headers], np.float32),
            "event_ID": np.array([header.event_number for header in chunk_generator.headers], np.uint16),
        }
        if chunk_generator.headers[0].event_weight > -1:
            header_level_info["event_weight"] = np.array([header.event_weight for header in chunk_generator.headers], np.float32)

        # Cross section info
        if chunk_generator.cross_section:
            # Even though this is a dataset level quantity, we need to match the structure in order to zip them together for storage.
            # Since we're repeating the value, hopefully this will be compressed effectively.
            header_level_info["cross_section"] = np.full_like(header_level_info["event_plane_angle"], chunk_generator.cross_section.value)
            header_level_info["cross_section_error"] = np.full_like(header_level_info["event_plane_angle"], chunk_generator.cross_section.error)

        # Assemble all of the information in a single awkward array and pass it on.
        yield ak.zip(
            {
                # Header level info
                **header_level_info,
                # Particle level info
                "particle_ID": ak.values_astype(array_with_events[:, :, 1], np.int32),
                # We're only considering final state hadrons or partons, so status codes are limited to a few values.
                # -1 are holes, while >= 0 are signal particles (includes both the jet signal and the recoils).
                # So we can't differentiate the recoil from the signal.
                "status": ak.values_astype(array_with_events[:, :, 2], np.int8),
                "E": ak.values_astype(array_with_events[:, :, 3], np.float32),
                "px": ak.values_astype(array_with_events[:, :, 4], np.float32),
                "py": ak.values_astype(array_with_events[:, :, 5], np.float32),
                "pz": ak.values_astype(array_with_events[:, :, 6], np.float32),
                # We could skip eta and phi since we can always recalculate them. However, since we've already parsed
                # them, we may as well pass them along.:w
                "eta": ak.values_astype(array_with_events[:, :, 7], np.float32),
                "phi": ak.values_astype(array_with_events[:, :, 8], np.float32),
            },
            depth_limit=1,
        )


def full_events_to_only_necessary_columns_E_px_py_pz(arrays: ak.Array) -> ak.Array:
    columns_to_drop = ["eta", "phi"]
    columns_to_keep = [field for field in ak.fields(arrays) if field not in columns_to_drop]
    return ak.zip(
        {
            column: arrays[column] for column in columns_to_keep
        }, depth_limit=1
    )


def parse_to_parquet(base_output_filename: Union[Path, str], store_only_necessary_columns: bool,
                     input_filename: Union[Path, str], events_per_chunk: int, parser: str = "pandas",
                     max_chunks: int = -1, compression: str = "zstd", compression_level: Optional[int] = None) -> None:
    """ Parse the JETSCAPE ASCII and convert it to parquet, (potentially) storing only the minimum necessary columns.

    Args:
        base_output_filename: Basic output filename. Should include the entire path.
        store_only_necessary_columns: If True, store only the necessary columns, rather than all of them.
        input_filename: Filename of the input JETSCAPE ASCII file.
        events_per_chunk: Number of events to be read per chunk.
        parser: Name of the parser. Default: "pandas".
        max_chunks: Maximum number of chunks to read. Default: -1.
        compression: Compression algorithm for parquet. Default: "zstd". Options include: ["snappy", "gzip", "ztsd"].
            "gzip" is slightly better for storage, but slower. See the compression tests and parquet docs for more.
        compression_level: Compression level for parquet. Default: `None`, which lets parquet choose the best value.
    Returns:
        None. The parsed events are stored in parquet files.
    """
    # Validation
    base_output_filename = Path(base_output_filename)
    # Setup the base output directory
    base_output_filename.parent.mkdir(parents=True, exist_ok=True)
    # We will check which fields actually exist when writing.
    possible_fields_containing_floats = ["event_plane_angle", "event_weight", "cross_section", "cross_section_error", "px", "py", "pz", "E"]

    for i, arrays in enumerate(read(filename=input_filename, events_per_chunk=events_per_chunk, parser=parser)):
        # Reduce to the minimum required data.
        if store_only_necessary_columns:
            arrays = full_events_to_only_necessary_columns_E_px_py_pz(arrays)
        else:
            # To match the steps taken when reducing the columns, we'll re-zip with the depth limited to 1.
            # As of April 2021, I'm not certainly this is truly required anymore, but it may be needed for
            # parquet writing to be successful (apparently parquet couldn't handle lists of structs sometime
            # in 2020. The status in April 2021 is unclear, but not worth digging into now).
            arrays = ak.zip(
                dict(zip(ak.fields(arrays), ak.unzip(arrays))),
                depth_limit = 1
            )

        # If converting in chunks, add an index to the output file so the chunks don't overwrite each other.
        if events_per_chunk > 0:
            suffix = base_output_filename.suffix
            output_filename = (base_output_filename.parent / f"{base_output_filename.stem}_{i:02}").with_suffix(suffix)
        else:
            output_filename = base_output_filename

        # Optimize the output
        # Additional parquet options are based on https://stackoverflow.com/a/66854439/12907985
        # byte_stream_fields apparently only work for float fields. Other fields should be handled
        # by use_dictionary. Apparently it can't handle this automatically, we so we have to define it
        # ourselves. This is a bit brittle if fields change, but they don't change so often, and
        # it's simpler than parsing field types, so it should be fine for now.
        byte_stream_fields = [field for field in ak.fields(arrays) if field in possible_fields_containing_floats]
        dict_fields = [field for field in ak.fields(arrays) if field not in possible_fields_containing_floats]
        # logger.debug(f"dict_fields: {dict_fields}")
        # logger.debug(f"byte_stream_fields: {byte_stream_fields}")

        # Parquet with zlib seems to do about the same as ascii tar.gz when we drop unneeded columns.
        # And it should load much faster!
        ak.to_parquet(
            arrays, output_filename,
            compression=compression, compression_level=compression_level,
            explode_records=False,
            # Additional parquet options are based on https://stackoverflow.com/a/66854439/12907985
            #use_dictionary=True,
            #use_byte_stream_split=True,
            use_dictionary=dict_fields,
            use_byte_stream_split=byte_stream_fields,
        )

        # Break now so we don't have to read the next chunk.
        if (i + 1) == max_chunks:
            break


if __name__ == "__main__":
    #read(filename="final_state_hadrons.dat", events_per_chunk=-1, base_output_filename="skim/jetscape.parquet")
    for pt_hat_range in ["7_9", "20_25", "50_55", "100_110", "250_260", "500_550", "900_1000"]:
        print(f"Processing pt hat range: {pt_hat_range}")
        directory_name = "OutputFile_Type5_qhatA10_B0_5020_PbPb_0-10_0.30_2.0_1"
        filename = f"JetscapeHadronListBin{pt_hat_range}"
        parse_to_parquet(
            base_output_filename=f"skim/{filename}.parquet",
            store_only_necessary_columns=True,
            input_filename=f"/alf/data/rehlers/jetscape/osiris/AAPaperData/MATTER_LBT_RunningAlphaS_Q2qhat/{directory_name}/{filename}_test.out",
            events_per_chunk=20,
            #max_chunks=3,
        )
