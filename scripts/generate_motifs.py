# script for generating motifs from the CNN

from scars import scars_queries
import keras
from keras.models import load_model
import keras_genomics
import pyranges as pr
import pybedtools
import pandas as pd
import numpy as np
from Bio import motifs
import sys

ancestry = sys.argv[1]
cnn_file = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/" + ancestry + "/cnn/cnn_25.h5"
motif_output_folder = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/" + ancestry + "/motifs/"

# load cnn and cut off all but first layer
cnn = load_model(cnn_file, custom_objects={'RevCompConv1D': keras_genomics.layers.RevCompConv1D})
n_layers = len(cnn.layers)
for _ in range(n_layers-1): 
	cnn.pop()

# load test set of sequences
training_loci = pr.read_bed("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/" + ancestry + "/training_sites/training_loci.bed", nrows=int(1e6))
n_seqs = training_loci.length

# query sequences
genome = pybedtools.genome_registry.hg38
reference_fasta = "/rds/project/who1000-1/rds-who1000-cbrc/ref/UCSC/hg38/hg38.fa"
flank = cnn.input_shape[1]//2
sequence = scars_queries.query_sequence(training_loci, flank, genome, reference_fasta)

# generate predictions for test sequences
signals = cnn.predict(sequence)
is_considered = np.max(signals, axis=1) > 0
max_positions = np.argmax(signals, axis=1)

kernel_size = cnn.layers[0].kernel_size[0]
n_motifs = cnn.layers[0].filters
max_signal_sequences = np.zeros(shape=(n_seqs, 2*n_motifs, kernel_size, 4), dtype=np.int8)

for i in range(2*flank+1):
	for j in range(2*n_motifs):
		seq_ix = (max_positions[:, j] == i)		# which sequences have i as the max signal position
		
		subseq_lower_ix = max(kernel_size//2 - i, 0)			# keep some zeros if max signal is near the edge
		subseq_upper_ix = min(kernel_size + 2*(flank-1) - i, kernel_size)
		
		max_signal_sequences[seq_ix, j, subseq_lower_ix:subseq_upper_ix, :] = sequence[seq_ix, max(0, i-kernel_size//2):min(2*flank+1, i+kernel_size//2+1), :]

max_signal_sequences[~is_considered] = np.zeros(shape=(kernel_size,4))

# summarise over all sequences
PFMs = np.sum(max_signal_sequences, axis=0)
nucleotides = ['A','C','G','T']
for motif_ix in range(n_motifs):
	nucleotide_counts = dict(zip(nucleotides, PFMs[motif_ix].T))
	motif = motifs.Motif(counts=nucleotide_counts)
	fname = motif_output_folder + str(motif_ix) + '.png'
	motif.weblogo(fname=fname)				# requires internet connection



