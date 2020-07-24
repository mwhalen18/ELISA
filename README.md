# ELISA

This code allows calculation from raw ELISA data

DATA Format:
Expected well layout:

         1     2     3   4   5   6   7   8   9   10  11  12
     A  NSB   NSB   B0  B0  U6  U6  U14 U14 U22 U22 U30 U30 
     B  STD1  STD1  CH  CH  U7  U7  U15 U15 U23 U23 U31 U31
     C  STD2  STD2  CL  CL  U8  U8  U16 U16 U24 U24 U32 U32
     D  STD3  STD3  U1  U1  U9  U9  U17 U17 U25 U25 U33 U33
     E  STD4  STD4  U2  U2  U10 U10 U18 U18 U26 U26 U34 U34
     F  STD5  STD5  U3  U3  U11 U11 U19 U19 U27 U27 U35 U35
     G  STD6  STD6  U4  U4  U12 U12 U20 U20 U28 U28 U36 U36
     H  STD7  STD7  U5  U5  U13 U13 U21 U21 U29 U29 U37 U37

#CH = Control HIGH
#CL = CONTROL LOW

Code can be imported as a data.frame with column headers '1' - '12'. There should be no row names and any extra rows should be deleted so the code can read just the 12x8 grid.

Output is a summary of standards, summary of samples, and a standard curve.

Its clunky and highly personalized for me and my needs. May optimize it in the future.
