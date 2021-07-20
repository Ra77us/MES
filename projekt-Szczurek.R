# Zadanie obliczeniowe - Michal Szczurek
install.packages("ggplot2")
library(ggplot2)
n = 4

# rysowanie czegokolwiek zeby potem zobaczyc czy na pewno wykres sie odswiezyl
plot(0,0)

mes(n)


##################

#calkowanie numeryczne

# funkcja przeksztalcajaca funkcje, by mozliwe bylo calkowanie na przdziale 
#innym niz (-1,1)
standarize <- function(f,x,a,b){
  return(f((b-a)/2*x+(a+b)/2)) 
}
#calkowanie numeryczne - 5 punktow kwadratury
integral <- function(f,a,b){
  
  x1 = 0
  x2 = -(1/3) * sqrt(5-2*sqrt(10/7))
  x3 = (1/3) * sqrt(5-2*sqrt(10/7))
  x4 = -(1/3) * sqrt(5 + 2*sqrt(10/7))
  x5 = (1/3) * sqrt(5 + 2*sqrt(10/7))
  return((b-a)/2 * (standarize(f,x1,a,b) * (128/225) +
                      standarize(f,x2,a,b) * ((322 +13*sqrt(70))/900) +
                      standarize(f,x3,a,b) * ((322 +13*sqrt(70))/900) +
                      standarize(f,x2,a,b) * ((322 - 13*sqrt(70))/900) +
                      standarize(f,x3,a,b) * ((322 - 13*sqrt(70))/900)))
  
  
  " Wersja dla 2 pkt
  x1 = 1 / sqrt(3)
  x2 = -1/sqrt(3)
  return((b-a)/2*(standarize(f,x1,a,b) + standarize(f,x2,a,b)))
  "
  
}

#####################
# MES
mes <- function(n){
  PI <- 3.14159
  G <- 6.674302 * 10^(-11)
  h <- 3/n
  
  # L(v)
  l_func <- function(v,i,h){
    const <- 4 * PI * G
    # ograniczam granice calkowania do przedzialu, na ktorym funkcja jest niezerowa
    # zwieksza to dokladnosc calkowania numerycznego
    a <- (i-1)*h
    b <- (i+1)*h
    a <- min(a,1)
    b <- max(b,2)
    if (a>b){
      return (0)
    }
    return (const *  integral(v,a,b))
  }
  
  
  # fukcja zwraca funkcje e_i(x)
  e <- function(i, h){
    
    x_i <- i * h
    
    func <- function(x){
      
      if (x_i-h < x && x <= x_i){
        return ((x - x_i + h)/h)
      }
      
      else if (x_i < x && x < x_i+h){
        return ((x_i - x + h)/h)
      }
      else{
        return (0)
      }
      
    }
    return(func)
  }
  
  # funckja zwraca pochaodna e_i(x())
  der_e <- function(h,i){
    x_i <- i * h
    
    func <- function(x){
      
      if (x_i-h < x && x <= x_i){
        return (1/h)
      }
      
      else if (x_i < x && x < x_i+h){
        return (-1/h)
      }
      else{
        return (0)
      }
      
    }
    return (func)
  }
  
  # Pochodna u z daszkiem
  # u z daszkiem to 5 * e_0 + 4 * e_n
  u_func <- function(h){
    
    func <- function(x){
      
      if (x >= 0 && x < h){
        return (-5/h)
      }
      
      else if (x> 3 - h && x < 3)
        return (4/h)
      
      else{
        return (0)
      }
    }
    
    return (func)
  }
  
  # macierz wartosci B(e_i,e_j)
  b_matrix <- matrix(nrow = n-1, ncol = n-1)
  
  # wektor wartosci l z daszkiem 
  l_vector <- vector(mode="numeric",length = n-1)
  
  # funkcja dodaje czlon B do L wynikajacy z shiftu,
  # tworzac ja wykorzystalem fakt, ze u z daszkiem
  # jest niezerowa tylko w [0,h] i [3-h,3] (wiec pochodna tez)
  shift_correction <- function(h,i,start){
    u <- u_func(h)  # pochodna u z daszkiem
    v <- der_e(h,i)
    func <- function(x){
      return (-(u(x) * v(x)))
    }
    # korzystam z faktu ze wartosc != 0 tylko w 2 przedzialach
    if (start){
      return (integral(func, 0, h))
    }
    if (!start){
      return (integral(func, 3 - h, 3))
    }
    
  }
  
  # wypelnienie wektora l z daszkiem 
  for (i in c(1:(n-1))){
    
    l_vector[i] <- l_func(e(i,h),i,h)
  }
  
  # uwzglednienie shiftu
  l_vector[1] <- l_vector[1] - shift_correction(h,1,TRUE)
  l_vector[n-1] <- l_vector[n-1] - shift_correction(h,n-1, FALSE)
  
  # funkcja obliczajaca wypelnienie macierzy w miejscu (i,j)
  fill_matrix <- function(i,j,h){
    
    if (abs(i-j)>= 2){
      return (0)
    }
    func <- function(x){
      return (der_e(h,i)(x) * der_e(h,j)(x))
    }
    if (i == j){
      return (-integral(func,(i-1)*h,(i+1)*h))
    }
    
    a <- min(i,j)
    b <- max(i,j)
    
    return (-integral(func,a*h, b*h))
  }
  
  # wypelnienie macierzy
  for (i in c(1:(n-1))){
    for (j in c(1:(n-1))){
      b_matrix [i,j] <- fill_matrix(i,j,h)
    }
  }
  
  # policzenie wektora u
  u_vector <- solve(b_matrix, l_vector)
  
  # wyznaczenie przykladowych punktow do wykresu
  y <- c()
  x <- c()
  
  num_of_pts <- 101
  for (i in c(0:num_of_pts)){
    
    x_val = i * 3 / num_of_pts
    y_val = 0 
    for (j in 1:length(u_vector)){
      y_val = y_val + u_vector[j] * (e(j,h)(x_val))
      
    }
    # uwzglednienie shiftu
    y_val = y_val + 5 * ((e(0,h))(x_val)) + 4 * ((e(n,h))(x_val)) 
    x = append(x, x_val)
    y = append(y, y_val)
  }
  
  res = data.frame(x, y)
  #print(b_matrix)
  ggplot(data=res, aes(x=x, y=y, group = 1)) +
    geom_line() 
   
}


