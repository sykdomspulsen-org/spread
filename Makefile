# h/t to @jimhester and @yihui for this parse block:
# https://github.com/yihui/knitr/blob/dc5ead7bcfc0ebd2789fe99c527c7d91afb3de4a/Makefile#L1-L4
# Note the portability change as suggested in the manual:
# https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Writing-portable-packages
export PKGNAME=`sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION`
export PKGVERS=`sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION`
export PKGTARBALL=$(PKGNAME)_$(PKGVERS).tar.gz
export DATETIME=`date +%Y-%m-%d\ %H:%M:%S`
export DATETIMEUTC=`date -u +%Y-%m-%d\ %H:%M:%S`
export DATE=`date +%Y.%-m.%-d`

all: check

fix_description_date:
	sed -i "s/^Version: .*/Version: $(DATE)/" DESCRIPTION
	sed -i '/Date\/Publication:/d' DESCRIPTION # delete if exists
	echo "Date/Publication: $(DATETIMEUTC) UTC" >> DESCRIPTION #append to bottom
	cat DESCRIPTION

install_deps:
	Rscript \
	-e 'if (!requireNamespace("remotes")) install.packages("remotes")' \
	-e 'remotes::install_deps(dependencies = TRUE, upgrade = "never")'

.ONESHELL:
build:
	Rscript \
		-e 'if (!requireNamespace("remotes")) install.packages("remotes")' \
		-e 'remotes::install_deps(dependencies = TRUE, upgrade = "never")'
	R CMD build .

.ONESHELL:
check:
	Rscript \
		-e 'if (!requireNamespace("remotes")) install.packages("remotes")' \
		-e 'remotes::install_deps(dependencies = TRUE, upgrade = "never")'
	R CMD check --no-manual $(PKGNAME)_$(PKGVERS).tar.gz

	cd *.Rcheck

	if grep -Fq "WARNING" 00check.log
	then
		# code if found
		exit 1
	else
		# code if not found
		echo "NO WARNINGs"
	fi

	if grep -Fq "ERROR" 00check.log
	then
		# code if found
		exit 1
	else
		# code if not found
		echo "NO ERRORs"
	fi

install: install_deps build
	R CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz

# this happens outside docker
.ONESHELL:
drat_update:
	cd /mnt/n/sykdomspulsen_config/drat
	git config user.name "Sykdomspulsen"
	git config user.email "sykdomspulsen@fhi.no"
	git config push.default simple
	git checkout gh-pages
	git pull

# this happens inside docker
.ONESHELL:
drat_insert:
	Rscript -e 'drat::insertPackage(fs::dir_ls("/rpkg/", regexp=".tar.gz$\"), repodir = "/drat")'

# this happens outside of docker
.ONESHELL:
drat_push:
	sed -i "/## News/a - **$(PKGNAME) $(PKGVERS)** (linux) inserted at $(DATETIME)" /mnt/n/sykdomspulsen_config/drat/README.md
	sed -i '1001,$ d' /mnt/n/sykdomspulsen_config/drat/README.md # only keep first 1000 lines of readme
	git -C /mnt/n/sykdomspulsen_config/drat add -A
	git -C /mnt/n/sykdomspulsen_config/drat commit -am "Jenkins $(PKGNAME) $(PKGVERS)" #Committing the changes
	git -C /mnt/n/sykdomspulsen_config/drat push -f origin gh-pages #pushes to master branch

.ONESHELL:
drat_prune_history:
	cd /tmp
	git clone "git@github.com:folkehelseinstituttet/drat.git"
	cd drat
	git config user.name "Sykdomspulsen"
	git config user.email "sykdomspulsen@fhi.no"
	git config push.default simple
	git checkout gh-pages

	git checkout --orphan latest_branch
	git add -A
	git commit -am "Cleaning history" #Committing the changes
	git branch -D gh-pages #Deleting master branch
	git branch -m gh-pages #renaming branch as master
	git -C /tmp/drat push -f origin gh-pages #pushes to master branch

# this happens inside of docker
.ONESHELL:
pkgdown_build:
	Rscript -e 'devtools::install("/rpkg", dependencies = TRUE, upgrade = FALSE); pkgdown::build_site("/rpkg")'

# this happens outside of docker:
pkgdown_deploy:
	git add .
	git commit -am "Pkgdown built"
	git subtree split --prefix docs -b gh-pages # create a local gh-pages branch containing the splitted output folder
	git push -f origin gh-pages:gh-pages # force the push of the gh-pages branch to the remote gh-pages branch at origin
	git branch -D gh-pages # delete the local gh-pages because you will need it: ref

clean:
	@rm -rf $(PKGNAME)_$(PKGVERS).tar.gz $(PKGNAME).Rcheck
