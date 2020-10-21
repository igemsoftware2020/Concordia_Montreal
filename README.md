# iGEM2020 - AstroBio
A web app that determines genes expression changes in microgravity

---------------------------------------------------------------------------------------------------------------------------------------


#### Mnagaing GitHub
  
  Cloning
  
 * Open terminal
 * Choose directory
 * clone the program ```git clone https://github.com/MaherHassanain/iGEM2020.git```
 
Process for Repository

1. Create new branch locally
  * Branch name is a condensed description of the task done
  * Use ```git checkout -b branch_name``` to create branch locally

2. Commit changes to branch
  * ```git add name_of_file_to_add``` to add file (* for all files)
  * ```git commit -m "commit_message"``` to commit changes

3. Push changes to branch
  * ```push remote origin branch_name```
  * ```git push origin```

4. Pull request
  * When the task is finished, open a pull request on the branch
  * Set a reviwer to review the code and merge 
  
---------------------------------------------------------------------------------------------------------------------------------------

#### Contributors

* Maher Hassanain - Team Leader, Full Stack Developer (Express js, Node js, JS, HTML, CSS) , Database Management (MongoDB), Software Integration and Launching (AWS)
* Benjamin Clark - Transcritpomics Analysis (R), Shiny App Developer 
* Hajar Elmouddene - Content creator, Frontend Developer (HTML/CSS)
* Grecia Olano - Content creator, Frontend Developer (HTML/CSS/JS)

---------------------------------------------------------------------------------------------------------------------------------------
 
#### OS Information

* Both the express app and the shiny app can run on vriaty of operating systems, including Windows, Linux, and macOS. 

---------------------------------------------------------------------------------------------------------------------------------------

#### Installation for AstroBio web App

Downlaod
* Nodejs https://nodejs.org/en/download/ version (v8 + is preferred)
    * Node version 8.9.4+ is preferred
    * NPM version 6.13.0 is preffered
* Mongodb https://www.mongodb.com/try 
    * mongodb version 4.2.1 is preffered

* Windows:
    * Ensure you have nom and node in terminal. You can check by writing the following commands in terminal
        * Node --version
        * npm -- version 
    * Create an empty folder, name it "db" in the following directory C:\data\db in your windows machine
    * Once the folder is created, go to terminal and navigate to the following direcotry C:\Program Files\MongoDB\Server\(version)\bin
    * Once in bin directory, write the command "mongod", mongodb is now hosted in default port 27017
    * In case you wish to run your mongodb in another port, please make sure to edit code and adjust port number accordengly 
    
* MacOS:
    * 
    
* Linux:

#### Installation for Shiny App 
