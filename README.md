# Lolipop_cycles_count
The GitHub repository contains a MATLAB implementation of T. R. Halford and K. M. Chugg's algorithm for counting cycles of girth, girth+2, and girth+4 in Low-Density Parity-Check (LDPC) codes.

The tool is provided as a MATLAB script and can be easily used to count the cycles of different girths in LDPC codes. The tool takes an LDPC code matrix as input and outputs the number of cycles of various girths in the code.

To use the tool, simply provide an LDPC code parity-check matrix as input to the MATLAB script. The output will provide the number of cycles of girth, girth+2, and girth+4 in the code. This information can be useful for designing LDPC codes with specific properties or for analyzing the performance of existing LDPC codes.

Example of use: 



H=qc2sparse('6_4_24.txt');


[g_ Ng_ Ng2_ Ng4_ Ng_per_u_ Ng2_per_u_ Ng4_per_u_ Ng_per_p_ Ng2_per_p_ Ng4_per_p_ ] = girth(H) 


g - girth


Ng - number of cycles of size girth


Ng2 - number of cycles of size girth+2


Ng4 - number of cycles of size girth+4


Ng_per_u-number of cycles per Variable nodes (VNs)


Ng_per_p-number of cycles per Check nodes (CNs)


girthstat count statistics for QC-LDPC codes


just put circulant size, number of VN and CN

T. R. Halford and K. M. Chugg, “An algorithm for counting short cycles in bipartite graphs,” IEEE Trans. Inform. Theory,
vol. 52, no. 1, 2006 and
https://www.slideshare.net/UsatyukVasiliy/enumerating-cycles-in-bipartite-graph-using-matrix-approach
