" vimrc

" let s:labs_files = substitute(system("find . | grep '\.tex$'"),"\n",' ','g')
" let s:bibs_files = substitute(system("find . | grep '\.bib$'"),"\n",' ','g')
" call texcomplete#SetFiles('labs',s:labs_files)
" call texcomplete#SetFiles('bibs',s:bibs_files)

call textoggle#Clear()
call textoggle#Set('beamer', 1)
call textoggle#Set('acr',    1)
call textoggle#Set('verb',   1)
call textoggle#Reload()

Spell 1
ConcealToggle 2
nnoremap <buffer> <leader>p :! open %:r.pdf<CR><CR>
nnoremap <buffer> <F5> :call tex#Compile("make",'jobstart')<CR>
nnoremap <buffer> <F6> :call tex#Compile("make",'termopen')<CR>
