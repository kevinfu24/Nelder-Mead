Nelder_Mead = function(x,f,tolerance,keep_x = FALSE){
  #browser()
  x_0 = t(x)
  X = t(optimsimplex::optimsimplex(method = "spendley", x0 = x_0)$newobj$x)
  
  if (length(X[1,])==1){
    X = t(X)
  }
  
  #Chek to see if dimensions fit
  if (length(X[1,]) != length(X[,1])+1){
    stop("Dimension is not n x n + 1")
  }
  
  #Dimension
  n = length(X[,1])
  
  if (keep_x == TRUE){
    x_list = list()
  }
  
  #Defining helper functions
  Sorter = function(X,f){
    n = length(X[,1])
    f_values = vector()
    for (i in 1:(n+1)){
      f_values = c(f_values,f(X[,i]))
    }
    Pre_Order = cbind(f_values, 1:(n+1))
    Ordering = Pre_Order[order(Pre_Order[,1]),][,2]
    X = X[,Ordering]
    if (n == 1){
      X = t(as.matrix(X))
    }
    return(list("X" = X, "f" = Pre_Order[Ordering,][,1]  , "order" = Ordering))
  }
  Centroid = function(X) {return((1/n) * rowSums(X[,1:n]))}
  if (n == 1){
    Centroid = function(X) {return(X[1])}
  }
  Reflection = function(X,x_bar,t) {return(x_bar + t*(X[,n+1] - x_bar))}
  
  Ordered = Sorter(X,f)
  x_bar = Centroid(Ordered$X)
  
  
  
  #Defining the sorted X, f(X) and x_bar, it 
  k = 0
  while(k < 1000 && ((norm(Ordered$X[,1]-Ordered$X[,n+1],type = "2") > tolerance) || abs(f(Ordered$X[,1]) - f(Ordered$X[,n+1])) > tolerance )){
    #Defining the sorted X, f(X) and x_bar, it 
    Ordered = Sorter(Ordered$X,f)
    x_bar = Centroid(Ordered$X)
    k = k + 1 
    
    if(keep_x == TRUE){
      x_list[[paste("x",k,sep='_')]] = Ordered$X
    }
    
    #defining x_{-1} and f(x_{-1})
    x_1 = Reflection(Ordered$X,x_bar,-1)
    f_1 = f(x_1)
    
    
    if (Ordered$f[1] <= f_1 &&  f_1 <= Ordered$f[n]){
      Ordered$f[n+1] = f_1
      Ordered$X[,n+1] = x_1 
      next 
    }
    else if (f_1 < Ordered$f[1]) {
      x_2 = Reflection(Ordered$X,x_bar,-2)
      f_2 = f(x_2)
      if (f_2 < f_1){
        Ordered$f[n+1] = f_2
        Ordered$X[,n+1] = x_2 
        next
      }else{
        Ordered$f[n+1] = f_1
        Ordered$X[,n+1] = x_1 
        next
      }
    }else if (f_1 >= Ordered$f[n]) {
      if(Ordered$f[n] <= f_1 && f_1 < Ordered$f[n+1]){
        x_half = Reflection(Ordered$X,x_bar,-0.5)
        f_half = f(x_half)
        if(f_half <= f_1){
          Ordered$f[n+1] = f_half
          Ordered$X[,n+1] = x_half 
          next
        }else{
          for (i in 2:(n+1)){
            Ordered$X[,i] = 0.5*(Ordered$X[,i]+Ordered$X[,1])
            Ordered$f[i] =  f(Ordered$X[,1])
            next
          }
        }
      }else if (f_1 >= Ordered$f[n+1]){
        x_half = Reflection(Ordered$X,x_bar,-0.5)
        f_half = f(x_half)
        if(f_half <= Ordered$f[n+1]){
          Ordered$f[n+1] = f_half
          Ordered$X[,n+1] = x_half
          next
        }else{
          for (i in 2:(n+1)){
            Ordered$X[,i] = 0.5*(Ordered$X[,i]+Ordered$X[,1])
            Ordered$f[i] =  f(Ordered$X[,1])
            next
          }
        }
      }
    }
  }
  if(keep_x == TRUE){
    return(list("X" = Ordered$X, "f" = Ordered$f, "Order" = Ordered$order, "X sequence" = x_list))
  }
  else{
    return(Ordered)
  }
}

#dir.create("examples")
#setwd("examples")
setwd("C:\\Users\\toke\\Documents\\examples")

#f = function(x) 100*(x[2] - x[1]^2)^2 + (1 - x[1])^2
f = function(x) 50*x[2]^2 + x[1]^2 + x[1]*x[2]


x_0 = c(-1.2, 1)
Result = Nelder_Mead(x_0,f,1e-10,keep_x = TRUE)

N = 100
xs = seq(-3, 2, length = N)
ys = seq(-3, 2, length = N)
zs = matrix(0, N, N)
for (i in 1:N) for (j in 1:N) zs[i, j] = f(c(xs[i], ys[j])) # slow
levels = c(.1, .3, 1:5, 10, 20, 30, 40, 50, 60, 80, 100, 500, 1000)
X_values = Result$`X sequence`
#X_values[[1]][,1]
#i=1
png(file="example%03d.png", width=800, height=800)
  for (i in 1:length(X_values)){
    contour(xs, ys, zs, levels = levels, col = "grey", title =title(main = paste("Iteration",i)))
    points(t(X_values[[i]][,1]),col= "blue")
    points(t(X_values[[i]][,2]), col = "blue")
    points(t(X_values[[i]][,3]), col = "blue")
    segments(t(X_values[[i]][,1])[1],t(X_values[[i]][,1])[2],t(X_values[[i]][,2])[1],t(X_values[[i]][,2])[2],col="blue")
    segments(t(X_values[[i]][,2])[1],t(X_values[[i]][,2])[2],t(X_values[[i]][,3])[1],t(X_values[[i]][,3])[2],col="blue")
    segments(t(X_values[[i]][,1])[1],t(X_values[[i]][,1])[2],t(X_values[[i]][,3])[1],t(X_values[[i]][,3])[2],col="blue")
  }
dev.off()

system('"C:\\Program Files\\ImageMagick-7.0.7-Q16\\convert.exe" -delay 15 -loop 0 *.png animation.gif')
file.remove(list.files(pattern=".png"))


