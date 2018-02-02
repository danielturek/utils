(setq backup-directory-alist `(("." . "~/.emacssave")))

;; this would allow aquamacs to detect rmd-mode automatically,
;; but the require commands cause an init-file error in bash emacs.
;;(setq load-path
;;    (append (list "/Users/dturek/utils/polymode/"
;;                  "/Users/dturek/utils/polymode/modes/")
;;        load-path))
;;(require 'poly-R)
;;(require 'poly-markdown)
;;(add-to-list 'auto-mode-alist '("\\.Rmd\\'" . poly-markdown+r-mode))
;;(add-to-list 'auto-mode-alist '("\\.rmd\\'" . poly-markdown+r-mode))

(defun rmd-mode ()
  "ESS Markdown mode for rmd files"
  (interactive)
  (setq load-path
    (append (list "/Users/dturek/github/utils/polymode/"
                  "/Users/dturek/github/utils/polymode/modes/")
        load-path))
  (require 'poly-R)
  (require 'poly-markdown)
  (poly-markdown+r-mode))

(defun rmarkdown-new-chunk ()
;;(defun rmarkdown-new-chunk (name)
  "Insert a new R chunk."
  (interactive)
  (insert "\n```{r }\n")
;;  (interactive "sChunk name: ")
;;  (insert "\n```{r " name "}\n")
  (save-excursion
    (newline)
    (insert "```\n")
    (previous-line)))

(defun rmarkdown-new-comment ()
  "Insert a new Rmarkdown comment."
  (interactive)
  (insert "\n<!--\n")
  (save-excursion
    (newline)
    (insert "-->\n")
    (previous-line)))

(defun rmarkdown-break-code-chunk ()
  "Insert a break into an Rmarkdown code chunk."
  (interactive)
  (insert "\n```\n")
  (save-excursion
    (newline)
    (insert "```{r }\n")
    (previous-line)))

(defun rmarkdown-weave-file ()
  "Run knitr on the current file and weave it as MD and HTML."
  (save-buffer)
  (interactive)
  (shell-command
   (format "/Users/dturek/github/utils/knit.sh -c %s"
       (shell-quote-argument (buffer-file-name)))))

(defun rmarkdown-tangle-file ()
  "Run knitr on the current file and tangle its R code."
  (interactive)
  (shell-command
   (format "/Users/dturek/github/utils/knit.sh -t %s"
       (shell-quote-argument (buffer-file-name)))))

(defun rmarkdown-commit ()  ;;(message)
  "Commit current buffer's repo"
  (interactive)
;;  (interactive "commit message: ")
  (shell-command
   ;(format "cd $(dirname %s); ~/github/utils/commit.sh "
   (format "cd $(dirname %s); FILEBASE=$(basename %s | cut -d. -f1);  git add $FILEBASE.*; git commit -m'.'; git push"
	   (shell-quote-argument (buffer-file-name))
	   (shell-quote-argument (buffer-file-name))
     )))

(defun rmarkdown-preview-file ()
  "Run knitr on the current file and display output in a browser."
  (save-buffer)
  (interactive)
  (shell-command
   (format "/Users/dturek/github/utils/knit.sh -b %s"
       (shell-quote-argument (buffer-file-name)))))

(global-set-key (kbd "C-c r") (quote R))
(global-set-key (kbd "C-c m") (quote rmd-mode))

(global-set-key (kbd "C-c n") (quote rmarkdown-new-chunk))
(global-set-key (kbd "C-c c") (quote rmarkdown-new-comment))
(global-set-key (kbd "C-c b") (quote rmarkdown-break-code-chunk))
(global-set-key (kbd "C-c w") (quote rmarkdown-weave-file))
(global-set-key (kbd "C-c p") (quote rmarkdown-preview-file))
(global-set-key (kbd "C-c t") (quote rmarkdown-tangle-file))
(global-set-key (kbd "C-c i") (quote rmarkdown-commit))


;; all this next stuff is for doing LaTeX docs from Aquamacs
;;(require 'tex-buf)  ;; causes a load-error when starting emacs from terminal
(defun TeX-command-default (name)
  "Next TeX command to use. Most of the code is stolen from `TeX-command-query'."
  (cond ((if (string-equal name TeX-region)
                 (TeX-check-files (concat name "." (TeX-output-extension))
                          (list name)
                          TeX-file-extensions)
               (TeX-save-document (TeX-master-file)))
             TeX-command-default)
            ((and (memq major-mode '(doctex-mode latex-mode))
                  (TeX-check-files (concat name ".bbl")
                           (mapcar 'car
                               (LaTeX-bibliography-list))
                           BibTeX-file-extensions))
             ;; We should check for bst files here as well.
             TeX-command-BibTeX)
            ((TeX-process-get-variable name
                           'TeX-command-next
                           TeX-command-Show))
            (TeX-command-Show)))


(defcustom TeX-texify-Show t "Start view-command at end of TeX-texify?" :type 'boolean :group 'TeX-command)
(defcustom TeX-texify-max-runs-same-command 5 "Maximal run number of the same command" :type 'integer :group 'TeX-command)

(defun TeX-texify-sentinel (&optional proc sentinel)
  "Non-interactive! Call the standard-sentinel of the current LaTeX-process.
If there is still something left do do start the next latex-command."
  (set-buffer (process-buffer proc))
  (funcall TeX-texify-sentinel proc sentinel)
  (let ((case-fold-search nil))
    (when (string-match "\\(finished\\|exited\\)" sentinel)
      (set-buffer TeX-command-buffer)
      (unless (plist-get TeX-error-report-switches (intern (TeX-master-file)))
    (TeX-texify)))))

(defun TeX-texify ()
  "Get everything done."
  (interactive)
  (let ((nextCmd (TeX-command-default (TeX-master-file)))
    proc)
    (if (and (null TeX-texify-Show)
         (equal nextCmd TeX-command-Show))
    (when  (called-interactively-p 'any)
      (message "TeX-texify: Nothing to be done."))
      (TeX-command nextCmd 'TeX-master-file)
      (when (or (called-interactively-p 'any)
        (null (boundp 'TeX-texify-count-same-command))
        (null (boundp 'TeX-texify-last-command))
        (null (equal nextCmd TeX-texify-last-command)))
    (mapc 'make-local-variable '(TeX-texify-sentinel TeX-texify-count-same-command TeX-texify-last-command))
    (setq TeX-texify-count-same-command 1))
      (if (>= TeX-texify-count-same-command TeX-texify-max-runs-same-command)
      (message "TeX-texify: Did %S already %d times. Don't want to do it anymore." TeX-texify-last-command TeX-texify-count-same-command)
    (setq TeX-texify-count-same-command (1+ TeX-texify-count-same-command))
    (setq TeX-texify-last-command nextCmd)
    (and (null (equal nextCmd TeX-command-Show))
         (setq proc (get-buffer-process (current-buffer)))
         (setq TeX-texify-sentinel (process-sentinel proc))
         (set-process-sentinel proc 'TeX-texify-sentinel))))))

(add-hook 'LaTeX-mode-hook '(lambda () (local-set-key (kbd "C-c C-c") 'TeX-texify))) 

