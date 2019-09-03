# Lolipop_cycles_count
T. R. Halford and K. M. Chugg Matlab implementation for LDPC codes cycles of girth, girth+2, girth+4 counting
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
