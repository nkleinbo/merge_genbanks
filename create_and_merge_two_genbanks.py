from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

def merge_genbank_files(file1, file2, output_file):
    # Step 1: Read the GenBank files
    record1 = SeqIO.read(file1, "genbank")
    record2 = SeqIO.read(file2, "genbank")
    
    # Step 2: Concatenate the sequences
    merged_sequence = record1.seq + record2.seq
    
    # Step 3: Adjust and combine the features
    merged_features = record1.features[:]
    offset = len(record1.seq)
    
    for feature in record2.features:
        # Adjust the feature location by adding the offset
        new_location = FeatureLocation(
            start=feature.location.start + offset,
            end=feature.location.end + offset,
            strand=feature.location.strand
        )
        # Create a new feature with the adjusted location
        new_feature = SeqFeature(
            location=new_location,
            type=feature.type,
            qualifiers=feature.qualifiers
        )
        merged_features.append(new_feature)
    
    # Step 4: Create a new SeqRecord with the concatenated sequence and combined features
    merged_record = SeqRecord(
        seq=merged_sequence,
        id=record1.id,
        name=record1.name,
        description=record1.description,
        features=merged_features,
        annotations=record1.annotations
    )
    
    # Step 5: Write the new SeqRecord to a GenBank file
    SeqIO.write(merged_record, output_file, "genbank")
    print(f"Merged GenBank file written to {output_file}")

#Create test objects / files and test it:

seq1 = Seq("ACCGCAGCGCAGTTACTAGCATCAGCATACATCGACTGATC")
seq2 = Seq("ACCGCAGCATGATGGTCATCAGCATACACCCCCCCCTCGACTGATC")

feat1 = SeqFeature(FeatureLocation(start = 5, end = 10), type="gene",  qualifiers={"gene": "example_gene"})
feat2 = SeqFeature(FeatureLocation(start = 8, end = 18), type="gene",  qualifiers={"gene": "example_gene"})

record1 = SeqRecord(
    seq=seq1,
    id="example_1",
    name="Example_1",
    description="An example SeqRecord number 1",
    features=[feat1, feat2],
    annotations = {"molecule_type": "DNA"}
)

record2 = SeqRecord(
    seq=seq2,
    id="example_2",
    name="Example_2",
    description="An example SeqRecord number 2",
    features=[feat1, feat2],
    annotations = {"molecule_type": "DNA"}
)

out1 = "example1.gbk"
out2 = "example2.gbk"

SeqIO.write(record1, out1, "genbank")
SeqIO.write(record2, out2, "genbank")

merge_genbank_files(out1, out2, "merged.gbk")

