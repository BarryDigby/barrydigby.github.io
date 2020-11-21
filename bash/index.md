---
title: Bash Basics
layout: page
permalink: /Bash/
---

Before beginning the tutorials and being given the keys to the car (access to Lugh) you are assumed to have a basic handle on bash. You should be able to comfortably navigate the terminal and perform basic commands. Bash is a fundamental skill in bioinformatics - one you will find yourself using on a daily basis. Best to get comfortable with it early.

In this section I will cover:

- [Navigating the file system](#navigating)
- [Organising files](#organising)
- [Deleting files](#deleting)

## Navigating the file system {#navigating}
Once you open a terminal, you will be brought to a directory that is predifined in your `.bashrc`. For most of you, this will be `~/` which is the same as `/home/username`. You can double check the directory you are in by typing `pwd`, short for `print working directory`.

To list the directories and files in your current directory listed by `pwd`, type `ls`. There are several useful variations of `ls` which I will point out to you in the example below.

To change to a different directory, use the command `cd`, short for `change directory`. Once inside a directory, if you want to go back to the parent directory, type `cd ../`. This command can be expanded upon, in the example below we enter the directory `/home/bdigby/perl/man`. To return to `/home/bdigby`, we go back two directories by typing `cd ../../`.

<center>
<img src="https://raw.githubusercontent.com/BarryDigby/BarryDigby.github.io/master/_gifs/basics_1.gif" width="100%" height="100%"/>
</center>

I used `ls -la` and `ls -ltr` in the terminal above. `ls -la` will print (from left to right) the user permissions of the file or directory, the owner of the file or directory, along with the size in bytes and time of creation. `ls -ltr` does the same thing, but it will order the results of `ls -la` so that the most recent files are at the bottom.

Bash is full of useful shortcuts, one of which is the wilcard denoted by and asterisk `*`. For example, if we wanted to list all text files in a directory with the extension `.txt`, we would use the command `ls *.txt` (or a variation of `ls` of your choice). Wildcard patterns are also referred to as `glob patterns`.

<center>
<img src="https://raw.githubusercontent.com/BarryDigby/BarryDigby.github.io/master/_gifs/wildcard.gif" width="100%" height="100%"/>
</center>

In the clip above I used the wildcard pattern to return only `.txt` files in the directory of interest. Another useful command to see the structure of directories in the current directory is the command `tree` as shown in the clip. This will return the contents of every directory (and its subdirectories) in the current directory.

## Organising Files {#organising}
Now that you can navigate and list the contents of directories, you need to be able to copy and move them. Bash accomplishes this with the `cp` and `mv` commands, which stand for `copy` and `move`, respectively. The commands follow this simple structure:

```bash
cp "path_to_file/file" "copy_to_destination"
mv "path_to_file/file" "move_to_destination"
```

The `cp` and `mv` commands also support glob patterns. To move all text files from the directory `/home/user/Download/Large_text_Download` to your Desktop, simply type `mv /home/user/Download/Large_text_Download/*.txt ~/Desktop`. Note that in order to copy and move files, you must have permissions to do so. This can be checked with `ls -la`. The image below gives a description of user permissions:

<center>
<img src="https://raw.githubusercontent.com/BarryDigby/BarryDigby.github.io/master/_images/file_permissions.png" width="100%" height="100%"/>
</center>

To the right of these codes will be the username of the file owner. You can ask the owner of the files to make them `rwx` for all users if it is appropriate to do so. Typically, large raw genomic data will be stored under private directories which require you to be added to a group in order to gain access (UKBB, TCGA, ICGC etc).

In order to move or copy entire directories, add the `-r` flag to the `cp` or `mv` command. For example: `cp -r ~/Downloads/python3 ~/usr/bin`


## Deleting files(#deleting)
To delete files, use the command `rm`. To delete directories, append the recursive flag `rm -r`. Bash will sometimes produce a message asking you if you are sure you want to delete a read/write protected file, prompting `(y/n)`, requiring the user to input an answer into the terminal. To bypass this message, you can use `rm -rf`, where the `-f` flag stands for `force`. The `-f` flag has the added benefit of attempting to remove the file - even if it does not exist - without producing a non 0 exit status.

A final note on `rm`. There is no 'Rubbish Bin'. Once it's gone, it's gone.

- [sed](#sed)
- [grep](#grep)
- [awk](#awk)

## Replacing strings with sed {#sed}
arandom words wronduinoi  wrkuciu erube eouihwcr wrckibw owouiwn owowcnroinwrc,,, wriubnc n crlnr;wrf974w wji wsknwhi8h wff4s8 okhs roi8 woh wd9oy8 9wchiuwhf9u i hwoih wcf49hg wc9ouhg 9w4gf 9wgh4 f9g 394gf 93g9dnbklj34e98yfd lkm3c;34rf...34 f083h 4f9n3 o94hf 3oef9 38y9fchpiegufrvarandom words wronduinoi  wrkuciu erube eouihwcr wrckibw owouiwn owowcnroinwrc,,, wriubnc n crlnr;wrf974w wji wsknwhi8h wff4s8 okhs roi8 woh wd9oy8 9wchiuwhf9u i hwoih wcf49hg wc9ouhg 9w4gf 9wgh4 f9g 394gf 93g9dnbklj34e98yfd lkm3c;34rf...34 f083h 4f9n3 o94hf 3oef9 38y9fchpiegufrv
arandom words wronduinoi  wrkuciu erube eouihwcr wrckibw owouiwn owowcnroinwrc,,, wriubnc n crlnr;wrf974w wji wsknwhi8h wff4s8 okhs roi8 woh wd9oy8 9wchiuwhf9u i hwoih wcf49hg arandom words wronduinoi  wrkuciu erube eouihwcr wrckibw owouiwn owowcnroinwrc,,, wriubnc n crlnr;wrf974w wji wsknwhi8h wff4s8 okhs roi8 woh wd9oy8 9wchiuwhf9u i hwoih wcf49hg wc9ouhg 9w4gf 9wgh4 f9g 394gf 93g9dnbklj34e98yfd lkm3c;34rf...34 f083h 4f9n3 o94hf 3oef9 38y9fchpiegufrvarandom words wronduinoi  wrkuciu erube eouihwcr wrckibw owouiwn owowcnroinwrc,,, wriubnc n crlnr;wrf974w wji wsknwhi8h wff4s8 okhs roi8 woh wd9oy8 9wchiuwhf9u i hwoih wcf49hg wc9ouhg 9w4gf 9wgh4 f9g 394gf 93g9dnbklj34e98yfd lkm3c;34rf...34 f083h 4f9n3 o94hf 3oef9 38y9fchpiegufrvarandom words wronduinoi  wrkuciu erube eouihwcr wrckibw owouiwn owowcnroinwrc,,, wriubnc n crlnr;wrf974w wji wsknwhi8h wff4s8 okhs roi8 woh wd9oy8 9wchiuwhf9u i hwoih wcf49hg wc9ouhg 9w4gf 9wgh4 f9g 394gf 93g9dnbklj34e98yfd lkm3c;34rf...34 f083h 4f9n3 o94hf 3oef9 38y9fchpiegufrvarandom words wronduinoi  wrkuciu erube eouihwcr wrckibw owouiwn owowcnroinwrc,,, wriubnc n crlnr;wrf974w wji wsknwhi8h wff4s8 okhs roi8 woh wd9oy8 9wchiuwhf9u i hwoih wcf49hg wc9ouhg 9w4gf 9wgh4 f9g 394gf 93g9dnbklj34e98yfd lkm3c;34rf...34 f083h 4f9n3 o94hf 3oef9 38y9fchpiegufrvarandom words wronduinoi  wrkuciu erube eouihwcr wrckibw owouiwn owowcnroinwrc,,, wriubnc n crlnr;wrf974w wji wsknwhi8h wff4s8 okhs roi8 woh wd9oy8 9wchiuwhf9u i hwoih wcf49hg wc9ouhg 9w4gf 9wgh4 f9g 394gf 93g9dnbklj34e98yfd lkm3c;34rf...34 f083h 4f9n3 o94hf 3oef9 38y9fchpiegufrvwc9ouhg 9w4gf 9wgh4 f9g 394gf 93g9dnbklj34e98yfd lkm3c;34rf...34 f083h 4f9n3 o94hf 3oef9 38y9fchpiegufrv




<br>

## Searching for strings with grep {#grep}
arandom words wronduinoi  wrkuciu erube eouihwcr wrckibw owouiwn owowcnroinwrc,,, wriubnc n crlnr;wrf974w wji wsknwhi8h wff4s8 okhs roi8 woh wd9oy8 9wchiuwhf9u i hwoih wcf49hg wc9ouhg 9w4gf 9wgh4 f9g 394gf 93g9dnbklj34e98yfd lkm3c;34rf...34 f083h 4f9n3 o94hf 3oef9 38y9fchpiegufrvarandom words wronduinoi  wrkuciu erube eouihwcr wrckibw owouiwn owowcnroinwrc,,, wriubnc n crlnr;wrf974w wji wsknwhi8h wff4s8 okhs roi8 woh wd9oy8 9wchiuwhf9u i hwoih wcf49hg wc9ouhg 9w4gf 9wgh4 f9g 394gf 93g9dnbklj34e98yfd lkm3c;34rf...34 f083h 4f9n3 o94hf 3oef9 38y9fchpiegufrvarandom words wronduinoi  wrkuciu erube eouihwcr wrckibw owouiwn owowcnroinwrc,,, wriubnc n crlnr;wrf974w wji wsknwhi8h wff4s8 okhs roi8 woh wd9oy8 9wchiuwhf9u i hwoih wcf49hg wc9ouhg 9w4gf 9wgh4 f9g 394gf 93g9dnbklj34e98yfd lkm3c;34rf...34 f083h 4f9n3 o94hf 3oef9 38y9fchpiegufrvarandom words wronduinoi  wrkuciu erube eouihwcr wrckibw owouiwn owowcnroinwrc,,, wriubnc n crlnr;wrf974w wji wsknwhi8h wff4s8 okhs roi8 woh wd9oy8 9wchiuwhf9u i hwoih wcf49hg wc9ouhg 9w4gf 9wgh4 f9g 394gf 93g9dnbklj34e98yfd lkm3c;34rf...34 f083h 4f9n3 o94hf 3oef9 38y9fchpiegufrvarandom words wronduinoi  wrkuciu erube eouihwcr wrckibw owouiwn owowcnroinwrc,,, wriubnc n crlnr;wrf974w wji wsknwhi8h wff4s8 okhs roi8 woh wd9oy8 9wchiuwhf9u i hwoih wcf49hg wc9ouhg 9w4gf 9wgh4 f9g 394gf 93g9dnbklj34e98yfd lkm3c;34rf...34 f083h 4f9n3 o94hf 3oef9 38y9fchpiegufrvarandom words wronduinoi  wrkuciu erube eouihwcr wrckibw owouiwn owowcnroinwrc,,, wriubnc n crlnr;wrf974w wji wsknwhi8h wff4s8 okhs roi8 woh wd9oy8 9wchiuwhf9u i hwoih wcf49hg wc9ouhg 9w4gf 9wgh4 f9g 394gf 93g9dnbklj34e98yfd lkm3c;34rf...34 f083h 4f9n3 o94hf 3oef9 38y9fchpiegufrvarandom words wronduinoi  wrkuciu erube eouihwcr wrckibw owouiwn owowcnroinwrc,,, wriubnc n crlnr;wrf974w wji wsknwhi8h wff4s8 okhs roi8 woh wd9oy8 9wchiuwhf9u i hwoih wcf49hg wc9ouhg 9w4gf 9wgh4 f9g 394gf 93g9dnbklj34e98yfd lkm3c;34rf...34 f083h 4f9n3 o94hf 3oef9 38y9fchpiegufrvarandom words wronduinoi  wrkuciu erube eouihwcr wrckibw owouiwn owowcnroinwrc,,, wriubnc n crlnr;wrf974w wji wsknwhi8h wff4s8 okhs roi8 woh wd9oy8 9wchiuwhf9u i hwoih wcf49hg wc9ouhg 9w4gf 9wgh4 f9g 394gf 93g9dnbklj34e98yfd lkm3c;34rf...34 f083h 4f9n3 o94hf 3oef9 38y9fchpiegufrv

## Column manipulation with awk {#awk}
arandom words wronduinoi  wrkuciu erube eouihwcr wrckibw owouiwn owowcnroinwrc,,, wriubnc n crlnr;wrf974w wji wsknwhi8h wff4s8 okhs roi8 woh wd9oy8 9wchiuwhf9u i hwoih wcf49hg wc9ouhg 9w4gf 9wgh4 f9g 394gf 93g9dnbklj34e98yfd lkm3c;34rf...34 f083h 4f9n3 o94hf 3oef9 38y9fchpiegufrvarandom words wronduinoi  wrkuciu erube eouihwcr wrckibw owouiwn owowcnroinwrc,,, wriubnc n crlnr;wrf974w wji wsknwhi8h wff4s8 okhs roi8 woh wd9oy8 9wchiuwhf9u i hwoih wcf49hg wc9ouhg 9w4gf 9wgh4 f9g 394gf 93g9dnbklj34e98yfd lkm3c;34rf...34 f083h 4f9n3 o94hf 3oef9 38y9fchpiegufrvarandom words wronduinoi  wrkuciu erube eouihwcr wrckibw owouiwn owowcnroinwrc,,, wriubnc n crlnr;wrf974w wji wsknwhi8h wff4s8 okhs roi8 woh wd9oy8 9wchiuwhf9u i hwoih wcf49hg wc9ouhg 9w4gf 9wgh4 f9g 394gf 93g9dnbklj34e98yfd lkm3c;34rf...34 f083h 4f9n3 o94hf 3oef9 38y9fchpiegufrvarandom words wronduinoi  wrkuciu erube eouihwcr wrckibw owouiwn owowcnroinwrc,,, wriubnc n crlnr;wrf974w wji wsknwhi8h wff4s8 okhs roi8 woh wd9oy8 9wchiuwhf9u i hwoih wcf49hg wc9ouhg 9w4gf 9wgh4 f9g 394gf 93g9dnbklj34e98yfd lkm3c;34rf...34 f083h 4f9n3 o94hf 3oef9 38y9fchpiegufrvarandom words wronduinoi  wrkuciu erube eouihwcr wrckibw owouiwn owowcnroinwrc,,, wriubnc n crlnr;wrf974w wji wsknwhi8h wff4s8 okhs roi8 woh wd9oy8 9wchiuwhf9u i hwoih wcf49hg wc9ouhg 9w4gf 9wgh4 f9g 394gf 93g9dnbklj34e98yfd lkm3c;34rf...34 f083h 4f9n3 o94hf 3oef9 38y9fchpiegufrvarandom words wronduinoi  wrkuciu erube eouihwcr wrckibw owouiwn owowcnroinwrc,,, wriubnc n crlnr;wrf974w wji wsknwhi8h wff4s8 okhs roi8 woh wd9oy8 9wchiuwhf9u i hwoih wcf49hg wc9ouhg 9w4gf 9wgh4 f9g 394gf 93g9dnbklj34e98yfd lkm3c;34rf...34 f083h 4f9n3 o94hf 3oef9 38y9fchpiegufrvarandom words wronduinoi  wrkuciu erube eouihwcr wrckibw owouiwn owowcnroinwrc,,, wriubnc n crlnr;wrf974w wji wsknwhi8h wff4s8 okhs roi8 woh wd9oy8 9wchiuwhf9u i hwoih wcf49hg wc9ouhg 9w4gf 9wgh4 f9g 394gf 93g9dnbklj34e98yfd lkm3c;34rf...34 f083h 4f9n3 o94hf 3oef9 38y9fchpiegufrvarandom words wronduinoi  wrkuciu erube eouihwcr wrckibw owouiwn owowcnroinwrc,,, wriubnc n crlnr;wrf974w wji wsknwhi8h wff4s8 okhs roi8 woh wd9oy8 9wchiuwhf9u i hwoih wcf49hg wc9ouhg 9w4gf 9wgh4 f9g 394gf 93g9dnbklj34e98yfd lkm3c;34rf...34 f083h 4f9n3 o94hf 3oef9 38y9fchpiegufrv
***
