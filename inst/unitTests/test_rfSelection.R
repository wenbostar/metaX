
test_rfSelection=function(){
    
    # case 1- users gives nothing to function
    checkException(rfSelection())
    
    # case 2- user gives only 1 argument instead of 2
    checkException(rfSelection(x=matrix(1:100,10,10)))
    
    # case 3- users gives uncorrect argument
    checkException(rfSelection(x=matrix(1:100,10,10),y=1:1000))
    
}
