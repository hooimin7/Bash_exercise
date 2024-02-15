ls -a ~ | grep .bash_profile # check if bash_profile exist
nano ~/.bash_profile # create bash_profile
# create a list of the numbers 1 through 5 called number.list
# all of these commands do the same thing
echo "1 2 3 4 5" | sed s/" "/"\n"/g > number.list # sed is a stream editor
echo "1 2 3 4 5" | tr " " "\n" > number.list # tr is a translate command
echo {1..5} | tr " " "\n" > number.list # brace expansion
seq 1 5 > number.list # seq is a sequence generator
for i in $(cat number.list)
do
echo "Welcome $i times"
done
##ctrl x to exit
##change permissions
# chmod +x ./mybash
chmod +x 
#run your script
sh mybash.sh