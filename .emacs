;; ____________________________________________________________________________
;; Aquamacs custom-file warning:
;; Warning: After loading this .emacs file, Aquamacs will also load
;; customizations from `custom-file' (customizations.el). Any settings there
;; will override those made here.
;; Consider moving your startup settings to the Preferences.el file, which
;; is loaded after `custom-file':
;; ~/Library/Preferences/Aquamacs Emacs/Preferences
;; _____________________________________________________________________________

;; Added by Package.el.  This must come before configurations of
;; installed packages.  Don't delete this line.  If you don't want it,
;; just comment it out by adding a semicolon to the start of the line.
;; You may delete these explanatory comments.
(package-initialize)

(setq backup-directory-alist `(("." . "~/.emacssave")))

;; this would allow aquamacs to detect rmd-mode automatically,
;; but the require commands cause an init-file error in bash emacs.
(setq load-path
    (append (list "~/github/utils/polymode/"
 		  "~/github/utils/polymode/modes/")
 	load-path))
(require 'poly-R)
(require 'poly-markdown)
(add-to-list 'auto-mode-alist '("\\.Rmd\\'" . poly-markdown+r-mode))
(add-to-list 'auto-mode-alist '("\\.rmd\\'" . poly-markdown+r-mode))
 
;;(unless (package-installed-p 'polymode)
;;  (package-install 'poly-markdown)
;;  (package-install 'poly-R))

(custom-set-variables
 '(one-buffer-one-frame-mode nil nil (aquamacs-frame-setup)))

(defun rmd-mode ()
  "ESS Markdown mode for rmd files"
  (interactive)
  (setq load-path
    (append (list "~/github/utils/polymode/"
                  "~/github/utils/polymode/modes/"
		  "~/github/utils/polymode/poly-R/")
        load-path))
  (require 'poly-R)
  (require 'poly-markdown)
  (poly-markdown+r-mode))

(defun rmarkdown-new-chunk ()
;;(defun rmarkdown-new-chunk (name)
  "Insert a new R chunk."
  (interactive)
  (insert "```{r }\n")
;;  (interactive "sChunk name: ")
;;  (insert "\n```{r " name "}\n")
  (save-excursion
    (newline)
    (insert "```\n")
    (previous-line)))

(defun rmarkdown-new-chunk-noeval ()
  "Insert a new R chunk with eval = FALSE."
  (interactive)
  (insert "```{r eval = FALSE}\n")
  (save-excursion
    (newline)
    (insert "```\n")
    (previous-line)))

(defun rmarkdown-new-comment ()
  "Insert a new Rmarkdown comment."
  (interactive)
  (insert "<!--\n")
  (save-excursion
    (newline)
    (insert "-->\n")
    (previous-line)))

(defun rmarkdown-break-code-chunk ()
  "Insert a break into an Rmarkdown code chunk."
  (interactive)
  (insert "```\n")
  (save-excursion
    (newline)
    (insert "```{r }\n")
    (previous-line)))

(defun rmarkdown-break-code-chunk-noeval ()
  "Insert a break into an Rmarkdown code chunk with eval = FALSE."
  (interactive)
  (insert "```\n")
  (save-excursion
    (newline)
    (insert "```{r eval = FALSE}\n")
    (previous-line)))

(defun rmarkdown-external-html-link ()
  "Insert template for external HTML link."
  (interactive)
  (insert "<a href=\"URL\" target=\"_blank\">TEXT</a>"))

(defun rmarkdown-external-html-link-blue ()
  "Insert template for external HTML link."
  (interactive)
  (insert "<a href=\"URL\" target=\"_blank\" style=\"color: blue\">TEXT</a>"))

(defun rmarkdown-insert-math-align ()
  "Insert align environment for math."
  (interactive)
  (insert "$$\\begin{align}\n")
  (save-excursion
    (newline)
    (insert "\\end{align}$$\n")
    (previous-line)))

(defun rmarkdown-insert-frac ()
  "Insert frac for math."
  (interactive)
  (insert "\\frac{}{}"))

(defun rmarkdown-weave-file ()
  "Run knitr on the current file and weave it as MD and HTML."
  (save-buffer)
  (interactive)
  (shell-command
   (format "~/github/utils/knit.sh -c %s"
       (shell-quote-argument (buffer-file-name)))))

(defun rmarkdown-tangle-file ()
  "Run knitr on the current file and tangle its R code."
  (interactive)
  (shell-command
   (format "~/github/utils/knit.sh -t %s"
       (shell-quote-argument (buffer-file-name)))))

(defun rmarkdown-commit ()  ;;(message)
  "Commit current buffer's repo"
  (interactive)
;;  (interactive "commit message: ")
  (shell-command
   ;(format "cd $(dirname %s); ~/github/utils/commit.sh "
   (format "cd $(dirname %s); FILEBASE=$(basename %s | cut -d. -f1);  git add $FILEBASE.*; git commit -m'updates'; git push"
	   (shell-quote-argument (buffer-file-name))
	   (shell-quote-argument (buffer-file-name))
     )))

(defun rmarkdown-preview-file ()
  "Run knitr on the current file and display output in a browser."
  (save-buffer)
  (interactive)
  (shell-command
   (format "~/github/utils/knit.sh -b %s"
       (shell-quote-argument (buffer-file-name)))))

(global-set-key (kbd "C-c r") (quote R))
(global-set-key (kbd "C-c m") (quote rmd-mode))

(global-set-key (kbd "C-c n") (quote rmarkdown-new-chunk))
(global-set-key (kbd "C-c N") (quote rmarkdown-new-chunk-noeval))
(global-set-key (kbd "C-c c") (quote rmarkdown-new-comment))
(global-set-key (kbd "C-c b") (quote rmarkdown-break-code-chunk))
(global-set-key (kbd "C-c B") (quote rmarkdown-break-code-chunk-noeval))
(global-set-key (kbd "C-c l") (quote rmarkdown-external-html-link))
(global-set-key (kbd "C-c L") (quote rmarkdown-external-html-link-blue))
(global-set-key (kbd "C-c a") (quote rmarkdown-insert-math-align))
(global-set-key (kbd "C-c f") (quote rmarkdown-insert-frac))
(global-set-key (kbd "C-c w") (quote rmarkdown-weave-file))
(global-set-key (kbd "C-c p") (quote rmarkdown-preview-file))
(global-set-key (kbd "C-c t") (quote rmarkdown-tangle-file))
(global-set-key (kbd "C-c i") (quote rmarkdown-commit))

;; next line is critically important!
;; makes it that in ESS, lines beginning with a single # comment, e.g.:
;; # this is a comment
;; are *not* indented 40 characters!!!  This is a *huge* improvement
(setq ess-fancy-comments nil)

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

(load-theme 'tango-dark)



