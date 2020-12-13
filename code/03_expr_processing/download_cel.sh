# while read f; 
# do  wget "$f"; 
# done < epithelium_list_cel_files.txt 

while IFS=',' read -r plat_type platform study fpath
do 
    mkdir -p "$plat_type/$platform/$study/";
    cd "$plat_type/$platform/$study/";
    wget "$fpath";
    cd ../../../;
done < "../list_f_to_download.csv"