for file in $2/*
do
	if [[ $file == *.json ]]
	then
		echo "processing $file..."
		{$1}/build/release/tools/simulator -i $3 -w -e $file --replicas 5 -o "$2"/"$2"_results/"${file##*/}"_quality.json --reload-binary
	fi
done

