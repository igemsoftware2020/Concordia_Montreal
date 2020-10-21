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

* Maher Hassanain - Team Leader, Full Stack Developer (Express js, Node js, AJAX, jQuery, HTML, CSS, Bootstrap) , Database Management (MongoDB), Software Integration and Launching (AWS)
* Benjamin Clark - Transcritpomics Analysis for data collection (R), Shiny App Developer 
* Hajar Elmouddene - Content creator, Frontend Developer (HTML/CSS/Bootstrap/jQuery)
* Grecia Olano - Content creator, Frontend Developer (HTML/CSS/Bootstrap/jQuery)

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
    * Mongodb can now be accessed via the mongo command
    
* MacOS:
    * Install Homebrew and XCode
    * You also have the option to install mongodb from terminal
         * brew update
        * brew install mongodb
        * mkdir -p /data/db
            * If permission was needed, 
                * sudo chown -R `id -un` /data/db
                * Then write password
         * mongod
         * Mongodb can now be accessed via the mongo command
    * Alternatevely, one can also install it as follow:
        * brew tap mongodb/brew
        * brew install mongodb-community@4.4 
        * brew services start mongodb-community@4.4
        * mongod --config /usr/local/etc/mongod.conf --fork
        * Mongodb can now be accessed via the mongo command
    
* Linux:
    * Create repository file that can be accessed with yum command after choosing which version to install
    * Assmuming we are going for version 4.4, create /etc/yum.d/mongodb-org-4.4.repo and add the following in the file:
    
      [mongodb-org-4.4]
      
      name=MongoDB Repository
      
      baseurl=https://repo.mongodb.org/yum/amazon/2013.03/mongodb-org/4.4/x86_64/
      
      gpgcheck=1
      
      enabled=1
      
      gpgkey=https://www.mongodb.org/static/pgp/server-4.4.asc
      
    * Then write the following commands
        * sudo yum install -y mongodb-org
        * sudo systemctl start mongod
        * You can verify status, sudo systemctl status mongod
        * Start mongodb service, sudo systemctl enable mongod
        * Mongodb can now be accessed via the mongo command
        * You can stop mongodb service, Stop >sudo systemctl stop mongod
        
---------------------------------------------------------------------------------------------------------------------------------------

#### Installation for Shiny App 

* Shiny app requires shiny server, shiny package, and many other packages in order to run successfully. 
* Therefore, we have created a setup_shiny.R which can be found in the repository in the following directory:
    * Concordia_Montreal/scripts/R/shinyApp/setup_shiny.R 
    * Running the setup_shiny.R code will install all dependent packages
* Now you can run app.R in Rstudio and explore the shiny app.

---------------------------------------------------------------------------------------------------------------------------------------

#### Installation for R scripts, transcriptomics analysis

* All R programs for data analytics are included in the following two directories: 
    * Concordia_Montreal/scripts/R/DEA/  
    * Concordia_Montreal/scripts/python
 * To set up all required packages for data analytics, run setup.R program which can be found in Concordia_Montreal/scripts/R/DEA/setup.R
