for file in $1/*
do
	if [[ $file == *.json ]]
	then
		echo "processing $file..."
		/global/homes/w/wadecap/ripples/build/release/tools/simulator -i $2 -w -e $file --replicas 5 -o "$1"/"$1"_results/"${file##*/}"_quality.json --reload-binary
	fi
done

