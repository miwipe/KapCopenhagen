for file in *.fa;do
bn=`basename $file .fa`
awk '/^>/ {gsub(/.fa(sta)?$/,"",FILENAME);printf(">%s\n",FILENAME);next;} {print}' $file > ${bn}_renamed.fa ; rm $file
done
