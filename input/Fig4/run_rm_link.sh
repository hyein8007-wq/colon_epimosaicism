cd /home/users/team_projects/hcolon_heterogeneity/01_LINE1/02_RNA/github/input/Fig4

for f in *; do 
	if [ -L "$f"  ] && [ -f "$f"  ]; then 
		cp -L "$f" "${f}.tmp" && mv "${f}.tmp" "$f"
	fi
done
