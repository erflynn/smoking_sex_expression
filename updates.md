
### Updat 9/8/2020
what I did:
- looked at a bit of the data for other tissues
  - most of it appears to be the same lab :/ 
- tried the sim data
  - this didnt show much of anything
- tried conormalizing with COCONUTS
  - this didnt appear to work as per looking at tSNE
- tried PVCA
  - also try tSNR
  - unclear if we can use JIVE for interactions
- tried coexpression networks with lioness
  - looks exciting! next steps:
    - use TF networks
    - assess
  - double check DE graph

### Update 7/15/2020
Added meta-analysis code as section 05a (cleaning), 05b (tables), 05c (exp), 05d (M-A).

Some next steps:

metadata:
- better study/sample mapping!, email PI
- double check metadata processing / lit comparison!!
- look into DGM IDs / repeated samples more
- fill in missing expression sex labels
- look at pm2.5
- possible supp analysis: look at age/genetic ancestry imputation

expression:
- where is the missing data? is it rescued elsewhere? 
- download from raw
- look into what is going on with duplicated data more!
 (specifically the data that is duplicated not by name)

visualization:
- tsne params
- try clustering
- consider pulling out cluster that may be tech artifact