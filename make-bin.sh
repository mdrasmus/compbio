



function install()
{
    dir1="$1"

    cd bin
    for file in "../$dir1"/*; do
	ln -s "$file" .
    done
    cd ..
}


# clean all links
for file in bin/*; do
    if (readlink "$file" > /dev/null); then
	rm "$file"
    fi
done

# install links
install bin-misc
install bin-phylogenomics
#install bin-compbio-misc

# hard code for Notung jar files
cd bin
ln -s ../bin-phylogenomics/Notung-2.6.jar .
ln -s ../bin-phylogenomics/Notung-2.7-Beta.jar .
cd ..
