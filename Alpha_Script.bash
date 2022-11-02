for data_file in ~/pythonProject1/Bird_MSAs/8782/4GI*
do
name=$(basename -s /home/will/pythonProject1/Bird_MSAs/8782/ $data_file)
#echo $data_file
output=$(basename -s .raw_alg.faa $name)
output=$output"_results.txt"
echo $output
python3 MET_program.py $data_file $output>>Bird_MSAs/processed/Combined
done
