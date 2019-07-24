for file in ../out/artifacts/*_jar/*.jar
do
	echo "Copying from ${file} to scripts/$(basename ${file})"
	cp ${file} scripts/$(basename ${file})
done
