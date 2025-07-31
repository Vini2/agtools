import gzip


class FastaParser:
    """
    A minimal, lightweight FASTA parser with on-demand sequence retrieval.

    This parser builds an index mapping sequence IDs to byte offsets in the file,
    allowing sequences to be fetched lazily without loading the entire FASTA into memory.
    Works with both plain-text FASTA and gzip-compressed FASTA (.gz).

    Attributes
    ----------
    file_path : str
        Path to the FASTA file (plain or gzipped).
    index : dict
        Mapping of sequence ID -> file offset for the header line.
    gzipped : bool
        True if the file is gzip-compressed.
    """

    def __init__(self, file_path, assembler="general", mapping=None):
        """
        Initialise the FastaParser and build an index for sequence IDs.

        Parameters
        ----------
        file_path : str
            Path to the FASTA file (.fasta or .fasta.gz).
        """
        self.file_path = file_path
        self.assembler = assembler
        self.mapping = mapping  # MEGAHIT
        self.index = {}
        self.gzipped = str(file_path).endswith(".gz")
        self._build_index()

    def _open(self, mode="rt"):
        """
        Open the FASTA file in text mode, supporting gzip if needed.

        Parameters
        ----------
        mode : str, optional
            File mode, default is 'rt' (read text).

        Returns
        -------
        file object
            An open file handle.
        """
        if self.gzipped:
            return gzip.open(self.file_path, mode)
        return open(self.file_path, mode)

    def _build_index(self):
        """
        Build an index mapping sequence IDs to byte offsets.

        For each header line starting with '>', store the current file position.
        This allows seeking to the start of a sequence later.
        """
        with self._open("rt") as f:
            pos = f.tell() if not self.gzipped else f.fileobj.tell()
            line = f.readline()
            while line:
                if line.startswith(">"):
                    seq_id = line[1:].strip().split()[0]
                    self.index[seq_id] = pos
                pos = f.tell() if not self.gzipped else f.fileobj.tell()
                line = f.readline()

    def get_sequence(self, seq_id):
        """
        Retrieve a DNA sequence by ID.

        Parameters
        ----------
        seq_id : str
            The sequence ID to fetch (matching the FASTA header without '>').

        Returns
        -------
        str or None
            The DNA sequence as a string, or None if the ID is not found.
        """
        seq_id = self.mapping[seq_id] if self.assembler == "megahit" else seq_id
        if seq_id not in self.index:
            return None
        seq_lines = []
        with self._open("rt") as f:
            if not self.gzipped:
                f.seek(self.index[seq_id])
            else:
                # For gzip, use fileobj.seek because f.seek is not fully random access
                f.fileobj.seek(self.index[seq_id])
            f.readline()  # skip header line
            for line in f:
                if line.startswith(">"):
                    break
                seq_lines.append(line.strip())
        return "".join(seq_lines)
