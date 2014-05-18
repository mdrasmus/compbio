



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
