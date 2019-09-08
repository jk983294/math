x <- 5
myenv <- new.env()
assign("x", "Homer", env = myenv) # assign in myenv
get("x", env = myenv) # get value from myenv
myenv$x <- "Hello" # assign in myenv
myenv$x # get value from myenv
ls()
ls(myenv)
x # global env's value not change

parent.env(myenv) # global env
