# awk to convert markdown file to C include for help.
{
    gsub(/\\/, "\\\\");  # Escape backslashes first
    gsub(/"/, "\\\"");   # Escape double quotes
    print "\"" $0 "\\n\""
}
