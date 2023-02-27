library(RColorBrewer)




# 
# sta <- list(type ="static",
#             dam_mortality =
#               list(
#                 e1 = 1,
#                 t1 = 5,
#                 e2 = 0,
#                 t2 = 45
#               ),
#             resource_rate =
#               list(
#                 r1 = 10,
#                 t1 = 5,
#                 r2 = 5,
#                 t2 = 45
#               )
#             ,
#             recruitment_max = 
#               list(
#                 rm1 = 0,
#                 t1 = 10,
#                 rm2 = 0,
#                 t2 = 40
#               )
# )

dyn <- list(type ="dynamic",
            dam_mortality =
              list(
                a = 1,
                k = 0.5),
            resource_rate =
              list(
                m = 10,
                u = -1.5,
                b = -0.05,
                start = 0.01,
                finish = 5),
            recruitment_max = list(
              l = 5,
              k = .5,
              x0 = 10
            ))


dm <- function(x)




# Mortality ---------------------------------------------------------------
dec_sig <- function(a,b,k,x) (-k/(1+exp(a+b*x)))+k
x <- seq(from = 1, to = 200, by = 1)
# b <- seq(-0.1, -1, by = -.1) #plot for paper
b <- seq(-3,-0.05,length.out =30) #scenario
k <- 1  #plot value is 1
a <- 0.5  # plot value is 0.5


res <- mapply(dec_sig, x = x, MoreArgs=list( b = b, k = k, a = a))

cols <- brewer.pal(n = 10, name = "Spectral")
par(mar = c(5, 5, 3, 8), xpd=TRUE)

matplot(x, t(res), col=cols, type="l", lty=1, lwd=2, xlab="Time", ylab="Dam Removal Mortality")
legend("right",  inset = c(- 0.25, 0), legend=b, title="Value of B", lwd=2, col=cols, bty="n")

# Resources ---------------------------------------------------------------

rebound_resource <- function(m,x,u,b) m-exp(-1*x)*(sin(u*x+b)/sin(b))
#test
m <- 10
u <- -1.5
x <- seq(0.01,10, length.out = 1000)
b <- seq(-.15, -0.05, length.out = 10) #plot
#b <- seq(-1, -0.05, length.out = 30) # scenario

res <- mapply(rebound_resource, x = x, MoreArgs=list( b = b, u = u,m = m))

cols <- brewer.pal(n = 10, name = "Spectral")
par(mar = c(5, 5, 3, 8), xpd=TRUE)

matplot(x, t(res), col=cols, type="l", lty=1, lwd=2, xlab="Time", ylab="Resource Rate")
legend("right",  inset = c(- 0.25, 0), legend=round(b,2), title="Value of b", lwd=2, col=cols, bty="n")

# Recruitment -------------------------------------------------------------

logis <- function(l,k,x,x0) ((l)/ (1+exp(-k*(x-x0))))

l = 0
k = 0.5
x = seq(1,50,1)
x0 = 10

Ls <- list(seq(1,5,length.out = 10)) #new asymptotes that range from 1 to 5

asym <- mapply(Ls, FUN = logis, k=k,x=x,x0=x0) #matrix of asymptotes over time

cols <- brewer.pal(n = 10, name = "Spectral")
par(mar = c(5, 5, 3, 8), xpd=TRUE)
matplot(x= t(matrix(rep(x,10),nrow=10,byrow = T)),y = t(asym),col=cols,type="l", lty=1, lwd=2, ylab = "Recruitment Max Scaling Amount",
        xlab = "Time")
legend("right",  inset = c(- 0.25, 0), legend=round(Ls[[1]],2), title="Value of l", lwd=2, col=cols, bty="n")
