## Create data

set.seed(56)
for(model_spec in c("asra1", "asra2")){
  for(condition in c("Asthma", "CHD")){
    for(sex in c("Persons", "Females", "Males")){
      
      N <- sample(100:1000, 1)
      
      alpha <- rnorm(1, 2,0.1)
      Beta <- rnorm(1, -3, 0.5)
      sigma <- runif(1, 0.01, 2)
      
      data <- data.frame(x = runif(N, -1,1),
                         sex = sex)
      data$y <- rnorm(N, alpha + Beta * data$x, sigma)
      

    }
    saveRDS(data, 
            paste0("C:/r_proj/QUTHPC_training/onHPC/data/", 
                   model_spec, '_', condition, ".rds"))
  }
}