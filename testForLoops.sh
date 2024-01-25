
for i in 2 2 2;
do
(
echo $i &
((j=i+1))
echo $j &
sleep 3 &
wait
)
done

