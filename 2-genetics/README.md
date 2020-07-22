# 2-genetics project
This is my R project for genetics
The graph results are in the folder ./deliverables/

## Steps and results
1) extract GSM matrix (.mtx) files as counts
2) reduced samples to just correspond to LCMV_GP66+ ("GSM3543444" & "GSM3543448") and TILs_CD44+_PD+ ("GSM3543446" "GSM3543449")
3) deleted any gene with below 1 CPM across all columns
4) TMM normalized gene epression and converted to Log-CPM
5) Used the treat method to find genes with significant cell type difference of LogFC greater than 2
  a) This left 12 genes with signficant difference across group
6) Used Independent 2 group Whitney-Man U test (a=.05) to test for differences between batches across each of the selected samples
  a) within groups, batches were all but twice (out of 12) significantly different from each other, suggesting significant batch differences.



### To-do: 
[x] add **/GSE124691/** to gitignore
[x] add 2-genetics/scratch/** to gitignore