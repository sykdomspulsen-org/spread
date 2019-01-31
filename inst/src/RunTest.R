con <- file("/tmp/computer","r")
COMPUTER_NAME <- readLines(con,n=1)
close(con)
Sys.setenv(COMPUTER=COMPUTER_NAME)

for(baseFolder in c("/data_clean","/results","/data_app")){
  files <- list.files(file.path(baseFolder,"noispiah"))
  if(length(files)>0){
    for(f in files) unlink(file.path(baseFolder,"noispiah",f))
  }
}

unlink(file.path("/junit","noispiah.xml"))
Sys.sleep(1)

a <- testthat:::JunitReporter$new()
a$start_reporter()
a$out <- file(file.path("/junit","noispiah.xml"), "w+")
a$start_context("noispiah")

output <- processx::run("Rscript","/r/noispiah/src/RunProcess.R", error_on_status=F, echo=T)
cat("\n\nstdout\n\n")
cat(output$stdout)
cat("\n\nstderr\n\n")
cat(output$stderr)

if(output$status==0){
  a$add_result("noispiah","RunAll",testthat::expectation("success","Pass"))
} else {
  a$add_result("noispiah","RunAll",testthat::expectation("error","Fail"))
}

a$end_context("noispiah")
a$end_reporter()
close(a$out)



