# Some differential equations in R

# Lorenz Model
library(deSolve)
library(rgl)

a <- -8/3; b <- -10; c <- 28
yini <- c(X = 1, Y = 1, Z = 1)
yini <- c(X = 2, Y = 2, Z = 2)

Lorenz <- function (t, y, parms) {
  with(as.list(y), {
    dX <- a * X+Y * Z
    dY <- b * (Y - Z)
    dZ <- -X * Y+c * Y-Z
    list(c(dX, dY, dZ)) })
}

times <- seq(from = 0, to = 100, by = 0.01)
out <- ode(y = yini, times = times, func = Lorenz,parms = NULL)
plot(out, lwd = 2)
plot(out[,"X"], out[,"Y"], type = "l", xlab = "X", ylab = "Y", main = "butterfly")

matplot(x = out[,"X"], y = out[,"Y"], type = "l", lwd = 2,
        lty = "solid", col = c("red", "blue", "black"),
        xlab = "time", ylab = "y", main = "rigidode")


plot3d(out[,"X"],out[,"Y"],out[,"Z"],col = "red",type = "l")

# Rigid Body Equations

yini <- c(1, 0, 0.9)
rigidode <- function(t, y, parms) {
  dy1 <- -2 * y[2] * y[3]
  dy2 <- 1.25* y[1] * y[3]
  dy3 <- -0.5* y[1] * y[2]
  list(c(dy1, dy2, dy3))
}

times <- seq(0,20, by = 0.01)
out <- ode (times = times, y = yini, func = rigidode,
            parms = NULL, method = rkMethod("rk45ck"))

head (out, n = 3)
matplot(x = out[,1], y = out[,-1], type = "l", lwd = 2,
        lty = "solid", col = c("red", "blue", "black"),
        xlab = "time", ylab = "y", main = "rigidode")
legend("bottomright", col = c("red", "blue", "black"),
       legend = c("y1", "y2", "y3"), lwd = 2)

plot3d(out[,-1])


# Arenstorf Orbits (Apollo Space Mission)

Arenstorf <- function(t, y, p) {
  D1 <- ((y[1] + mu1)^2 + y[2]^2)^(3/2)
  D2 <- ((y[1] - mu2)^2 + y[2]^2)^(3/2)
  dy1 <- y[3]
  dy2 <- y[4]
  dy3 <- y[1] + 2*y[4] - mu2*(y[1]+mu1)/D1 - mu1*(y[1]-mu2)/D2
  dy4 <- y[2] - 2*y[3] - mu2*y[2]/D1 - mu1*y[2]/D2
  return(list( c(dy1, dy2, dy3, dy4) ))
}
mu1 <- 0.012277471
mu2 <- 1 - mu1
yini <- c(y1 = 0.994, y2 = 0, dy1 = 0, dy2 = -2.00158510637908252240537862224)

times <- seq(from = 0, to = 18, by = 0.01)

out <- ode(func = Arenstorf, y = yini, times = times,parms = 0, method = "ode45")

# different initial conditions
yini2 <- c(y1 = 0.994, y2 = 0, dy1 = 0, dy2 = -2.0317326295573368357302057924)
out2 <- ode(func = Arenstorf, y = yini2, times = times,parms = 0, method = "ode45")

yini3 <- c(y1 = 1.2, y2 = 0,dy1 = 0, dy2 = -1.049357510)
out3 <- ode(func = Arenstorf, y = yini3, times = times,parms = 0, method = "ode45")

plot(out, out2, out3, which = c("y1", "y2"),mfrow = c(2, 2), col = "black", lwd = 2)

plot(out[ ,c("y1", "y2")], type = "l", lwd = 2,xlab = "y1", ylab = "y2", main = "solutions 1,2")
lines(out2[ ,c("y1", "y2")], lwd = 2, lty = 2)

plot(out3[ ,c("y1", "y2")], type = "l", lwd = 2, lty = 3,xlab = "y1", ylab = "y2", main = "solution 3")

plot(out[,"y1"],out[,"y2"])
plot(out2[,"y1"],out2[,"y2"])
plot(out3[,"y1"],out3[,"y2"])

plot(out[,"dy1"],out[,"dy2"])
plot(out2[,"dy1"],out2[,"dy2"])
plot(out3[,"dy1"],out3[,"dy2"])


# Pleiade 7 star system

pleiade <- function (t, Y, pars) {
  x <- Y[1:7]
  y <- Y[8:14]
  u <- Y[15:21]
  v <- Y[22:28]
  distx <- outer(x, x, FUN = function(x, y) x - y)
  disty <- outer(y, y, FUN = function(x, y) x - y)
  rij3 <- (distx^2 + disty^2)^(3/2)
  fx <- starMass * distx / rij3
  fy <- starMass * disty / rij3
  list(c(dx = u,
         dy = v,
         du = colSums(fx, na.rm = TRUE),
         dv = colSums(fy, na.rm = TRUE)))
  
}
starMass <- 1:7
yini<- c(x1= 3, x2= 3, x3=-1, x4=-3, x5= 2, x6=-2, x7= 2,
         y1= 3, y2=-3, y3= 2, y4= 0, y5= 0, y6=-4, y7= 4,
         u1= 0, u2= 0, u3= 0, u4= 0, u5= 0, u6=1.75, u7=-1.5,
         v1= 0, v2= 0, v3= 0, v4=-1.25, v5= 1, v6= 0, v7= 0)
print(system.time(
  out <- ode(func = pleiade, parms = NULL, y = yini,
             method = "adams", times = seq(0, 3, 0.01))))

par(mfrow = c(3, 3), mar = c(4, 4, 3, 2))
for (i in 1:7) {
  plot(out[,i+1], out[,i+8], type = "l",
       main = paste("star ",i), xlab = "x", ylab = "y")
  points (yini[i], yini[i+7])
}

#
plot(out[, 2:8], out[, 9:15], type = "p", cex = 0.5,main = "ALL", xlab = "x", ylab = "y")
text(yini[1:7], yini[8:14], 1:7)

plot(out[, c(2,8)], out[, c(9,15)], type = "p", cex = 0.5,main = "ALL", xlab = "x", ylab = "y")


#
matplot(out[,"time"], out[, c("u1", "u7")], type = "l",
        lwd = 2, col = c("black", "grey"), lty = 1,
        xlab = "time", ylab = "velocity", main = "stars 1, 7")
abline(v = c(1.23, 1.68), lty = 2)
legend("bottomright", col = c("black", "grey"), lwd = 2,
       legend = c("u1", "u7"))


# bouncing ball

yini <- c(height = 0, velocity = 10)
ball <- function(t, y, parms) {
  dy1 <- y[2]
  dy2 <- -9.8
  list(c(dy1, dy2))
}

rootfunc <- function(t, y, parms) y[1]
eventfunc <- function(t, y, parms) {
  y[1] <- 0
  y[2] <- -0.9*y[2]
  return(y)
}

times <- seq(from = 0, to = 20, by = 0.01)
out <- ode(times = times, y = yini, func = ball,
           parms = NULL, rootfun = rootfunc,
           events = list(func = eventfunc, root = TRUE))

plot(out, which = "height", lwd = 2,main = "bouncing ball", ylab = "height")
