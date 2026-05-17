from Bio import SeqIO
from stringzilla import Str, File
import gzip
import json


class FastaHolder:
    def __init__(self, filepath):
        self.filepath = str(
            filepath
        )  # path are given as a Path object so they need to be transformed to string
        self.sequences = list(SeqIO.parse(self.filepath, "fasta"))

    def get_sequence_lengths(self):
        return {record.id: len(record.seq) for record in self.sequences}

    def get_filepath(self):
        return self.filepath

    def get_dict_fasta(self):
        fastx_Str = fastx_io(self.filepath)
        return fasta_to_dict(fastx_Str)

    def save_fasta(self, output_path):
        with open(output_path, "w") as output_handle:
            for record_id, sequence in self.sequences.items():
                SeqIO.write(
                    SeqIO.SeqRecord(sequence, id=record_id), output_handle, "fasta"
                )


class FastqHolder:
    def __init__(self, filepath):
        self.filepath = str(
            filepath
        )  # path are given as a Path object so they need to be transformed to string
        self.sequences = parse_fastq_fasta(self.filepath)

    def get_sequence_lengths(self):
        return {record.id: len(record.seq) for record in self.sequences}

    def get_filepath(self):
        return self.filepath

    def get_dict_fasta(self):
        return {record.id: str(record.seq) for record in self.sequences}


# TODO put this in utils
def save_list_nodes_to_file(list_nodes, filename):
    with open(filename, "w") as file:
        json.dump(list_nodes, file, indent=4)


def fastx_io(path) -> Str:
    """
    Takes in the path of the fastx and returns the Str (stringzilla)
    to process with fasta_to_dict
    """
    # verify that the path is a path lib or a tst
    if isinstance(path, str):
        path = str(path)
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as f:
            return Str(f.read())
    else:
        return Str(File(path))


def fasta_to_dict(fastx: Str) -> dict:
    """
    Function to obtain the sequences from a fasta/(fastq ) file
    this functions uses stringzilla
    I could consider using the function from mappy TBD
    """
    # verify that it is a fastq
    fastx_lines = fastx.split(separator="\n", keepseparator=False)
    # length_file = len(fastx_lines) - 1
    if fastx[0] == "@":
        seqs = fastx_lines[1::4]  # seq are ever 4th line from the second
        ids = fastx_lines[::4]  #
        dict_fastq: dict[str, str] = {
            str(header[1:]): str(seq) for (header, seq) in zip(ids, seqs)
        }
        return dict_fastq
    else:
        # TODO verify that seq are written on 1 line e.g every second line should ahve a >header
        dict_fasta = {}
        id = ""
        try:
            for line in fastx_lines:
                if line != "":
                    if line[0] == ">":
                        if id != "":
                            dict_fasta[id] = "".join(dict_fasta.get(id))
                        id = str(line[1:])
                        dict_fasta[id] = []
                    else:
                        dict_fasta[id].append(str(line))
        finally:
            dict_fasta[id] = "".join(dict_fasta.get(id))
        # seqs = fastx_lines[1::2]
        # ids = fastx_lines[::2]
        # if ids[3,4,5] is not ">":

        # print(seqs)
        # dict_fasta   = {str(header[1:]): str(seq) for (header, seq) in zip(ids, seqs)}
        return dict_fasta

    # TODO check whether the last element of the seqs and id is not empty char


def parse_fastq_fasta(file_path):
    records = []
    # Determine if the file is compressed
    if file_path.endswith(".gz"):
        with gzip.open(file_path, "rt") as f:
            for record in SeqIO.parse(
                f, "fastq" if file_path.endswith(".fastq.gz") else "fasta"
            ):
                records.append(record)
    else:
        with open(file_path) as f:
            for record in SeqIO.parse(
                f, "fastq" if file_path.endswith(".fastq") else "fasta"
            ):
                records.append(record)

    return records


### legacy
def fasta_to_fastq(fasta_path: str, fastq_path: str):
    """
    Converts a FASTA file to a FASTQ file by adding default quality scores.

    Args:
        fasta_path (str): Path to the input FASTA file.
        fastq_path (str): Path to the output FASTQ file.
    """
    # Default quality score
    default_quality = "I" * 100  # Assuming sequence length is less than or equal to 100

    try:
        with open(fasta_path, "r") as fasta_file, open(fastq_path, "w") as fastq_file:
            for line in fasta_file:
                if line.startswith(">"):
                    # Write the header
                    header = line.strip()
                    fastq_file.write(header + "\n")
                else:
                    # Write the sequence
                    sequence = line.strip()
                    fastq_file.write(sequence + "\n")
                    # Write the quality scores
                    fastq_file.write("+" + "\n")
                    fastq_file.write(default_quality[: len(sequence)] + "\n")

        print(f"Conversion successful: {fasta_path} -> {fastq_path}")
    except FileNotFoundError:
        print(f"Error: File not found - {fasta_path}")
    except Exception as e:
        print(f"An error occurred during conversion: {e}")
