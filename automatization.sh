mkdir outputs

cd inputs

ls *.txt | sed 's/.txt//' > ../input

cd ..

parallel --jobs 1 "cp inputs/{}.txt enrichr_input.txt && python enrichr_api.py && mv annotationenrichr.tsv outputs/{}_output.tsv && rm enrichr_input.txt" :::: input
