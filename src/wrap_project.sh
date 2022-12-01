# Script to wrap and save project after it has been calculated
# also to make figures

# Create new dir in specified folder

mkdir -p $1
echo Project folder: $1

cp -R ../input $1
echo ../input folder copied to $1
rm ../input/secs_doys_sources_stations.txt
if test -d "$1/output"; then
    echo "$1/output exists, stopping..."
else
    if test -d ../output; then
        mv -n ../output $1
        echo "../output folder moved to $1"
    else
        echo "ERROR: file ../output doesn't exist"
    fi
fi