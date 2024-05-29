" vimrc

let s:labs_files = substitute(glob("./*.tex",), "\n", " ", "g")
let s:bibs_files = 'references.bib'
call texcomplete#SetFiles('labs', s:labs_files)
call texcomplete#SetFiles('bibs', s:bibs_files)

call textoggle#Clear()
call textoggle#Set('acr',  1)
call textoggle#Set('verb', 1)
call textoggle#Set('alg',  1)
call textoggle#Reload()

Spell 1
ConcealToggle 2
