pri=function()
{
  c=1:100
  for (i in 1:100)
  {
    if (i/3== round(i/3))
    {c[i]="Fizz"}
    if (i/5== round(i/5))
    {c[i]="Buzz"}
    if (i/15== round(i/15))
    {c[i]="FizzBuzz"}
  }
  return(c)
}