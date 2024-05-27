Compile　
Run "make" in the HRCM source code directory, the executable file hrcm is in the current directory.
Run "chmod 777 7za" to change the access permission of the 7-zip executables.

Usage
hrcm {compress | decompress}  -r {ref-file-path}{ [-t] {tar-file-path}|[-f] {filename} [percent]}
     {compress | decompress} is mode,  choose one of them according to requirement, required
     -r is the reference, the {ref-file-path} followed, required
     -t is the target, a single to-be-compressed file path {tar-file-path} followed, optional
     -f is the alternative option of -t, a set of to-be-compressed file paths included in {filename}, optional
     [percent] is the percentage of the second-level matching, default is 10, means 10% of sequences will be used for the second-level matching, optional when -f, illegal when -t

Output
1 compressed file named filename.7z 
2 decompressed file named filename.fasta

Example

    compress and decompress hg18_chr22.fa using hg17_chr22.fa as reference
    ./hrcm compress -r hg17_chr22.fa -t hg18_chr22.fa
    output: hg18_chr22.7z
	
    ./hrcm decompress -r hg17_chr22.fa -t hg18_chr22.7z
    output: hg18_chr22.fasta