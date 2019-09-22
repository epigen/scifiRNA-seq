
mkdir -p /scratch/users/$USER
cd /scratch/users/$USER
mkdir -p custom_demux/scRNA/
cd custom_demux/scRNA/
mkdir -p BSF_0666_HCW2YDRXX_1_samples BSF_0666_HCW2YDRXX_2_samples BSF_0666_HCW2YDRXX_3_samples BSF_0666_HCW2YDRXX_4_samples
cd BSF_0666_HCW2YDRXX_1_samples
ln -s /scratch/users/dbarreca/private/custom_demux/scRNA/BSF_0666_HCW2YDRXX/BSF_0666_HCW2YDRXX_1_samples/BSF_0666_HCW2YDRXX_1#*.bam .
cd ..
cd BSF_0666_HCW2YDRXX_2_samples
ln -s /scratch/users/dbarreca/private/custom_demux/scRNA/BSF_0666_HCW2YDRXX/BSF_0666_HCW2YDRXX_2_samples/BSF_0666_HCW2YDRXX_2#*.bam .
cd ..
cd BSF_0666_HCW2YDRXX_3_samples
ln -s /scratch/users/dbarreca/private/custom_demux/scRNA/BSF_0670_HFCK2DRXX/BSF_0670_HFCK2DRXX_1_samples/BSF_0670_HFCK2DRXX_1#*.bam .
rename BSF_0670_HFCK2DRXX_1# BSF_0666_HCW2YDRXX_3# *.bam
cd ..
cd BSF_0666_HCW2YDRXX_4_samples
ln -s /scratch/users/dbarreca/private/custom_demux/scRNA/BSF_0670_HFCK2DRXX/BSF_0670_HFCK2DRXX_2_samples/BSF_0670_HFCK2DRXX_2#*.bam .
rename BSF_0670_HFCK2DRXX_2# BSF_0666_HCW2YDRXX_4# *.bam
cd ..

