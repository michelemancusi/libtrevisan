# Trevisan-extractor
Implementation of Luca Trevisan's extractor

In this folder there are differents files: 

extractor.c is my code of the Trevisan Extractor; 

trevisan.c is the library in which there are the important functions to be able to execute the code (I thank James S. Plank - plank@cs.utk.edu - http://www.cs.utk.edu/~plank from whom I took the necessary functions to be able to do the calculations in small finite fields, while the functions to be able to do the calculations in finite fields of arbitrary size is my work);

trevisan.h is the header file;

polinomi_irriducibili.txt is the file in which the irreducible polynomials are saved in the form of strings of zeros and ones up to about the irreducible polynomial of GF (2 ^ 100). In order to make extractions and then calculations that require larger finite fields, just add the strings of zeros and ones corresponding to the irreducible polynomials you need at the end of this file;

source.txt and seed.txt are two sample files because the code requests in input a source file and a seed file. Obviously these two files will have to be replaced with your own source and seed files (always in the .txt extension containing zeros and ones).

When the extractor.c code is compiled and executed, it requires input data: the name of the source file, the name of the seed file, the value of the min-entropy per bit of the source and finally the desired error per bit (the smaller the error the more uniform the output string will be, but the longer the execution time will be).



P.S. If you have a segmentation fault after the input of the data, try to change sizeof(bool) with sizeof(int) in lines 191 and 197 of file extractor.c 
