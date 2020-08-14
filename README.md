# Tajima's D

## What is Tajima's D

 Tajima's D is the difference of two population genetics measures, the average number of pairwise differences and the total number of 
 segregating sites (bases where these differences occur). This is scaled to be behaving such that this would be the same as in a neutrally 
 evolving population of constant size.
 
The equation is as follows:

![alt text](https://s0.wp.com/latex.php?latex=D%3D%5Cfrac%7B%5Chat%7B%5Ctheta%7D_T+-+%5Chat%7B%5Ctheta%7D_W%7D%7B%5Csqrt%7B%5Chat%7BV%7D%28%5Chat%7B%5Ctheta%7D_T+-+%5Chat%7B%5Ctheta%7D_W%29%7D%7D&bg=ffffff&fg=111111&s=3)
 
 where 0t is the number of pairwise differences and 0w is the number of segregating sites.
 
## installation
```python setup.py install```

or
```pip install tajima```

## Running the script
```tajima -F fasta_file```

and to perform alignment using MUSCLE (if the FASTA is not of equal length):

```tajima -F fasta_file -M```
