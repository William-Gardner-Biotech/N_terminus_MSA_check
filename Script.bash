for data_file in ~/pythonProject1/mammalia/40674/*
do
name=$(basename -s /home/will/pythonProject1/mammalia/40674/ $data_file)
#echo $data_file
output=$(basename -s .raw_alg.faa $name)
output=$output"_results.txt"
#echo $output
python3 MET_program.py $data_file $output>mammalia/processed/$output
done
