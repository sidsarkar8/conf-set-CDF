### combing batch files Priyamvada gave me
### centred (x-mu)^2
### sample size = 500
### n_sim = 500

n2_seq = c(20,50,100,250,500)

for(n2 in n2_seq)
{
  ans = NULL
  for( i in 1:5){
    name_qq = paste("/Users/sidsarkar/Documents/Projects/conf_tail/other_methods/hydra/x2/",
                    n2,"_batch", i,".csv",sep = "")
    ans = rbind(ans, read.csv(name_qq))
  }
  
  name_ans = paste("/Users/sidsarkar/Documents/Projects/conf_tail/other_methods/hydra/x2/x2_", 
                   n2, ".csv", sep = "")
  write.csv(ans, name_ans)
}
