var = input("Do you want to clean up this directory? (y/n): ")

if var=="y":
    print( "Continued with code")
elif var=="n":
    print("Exiting script...")
    exit()
else:
    print("Invalid input - must be y/n")
    exit()
    
print("Code continued correctly")
