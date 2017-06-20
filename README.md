Installtion for softwares for tanglab NGS docker.

usage:

```bash
# avoiding Permission issue
mkdir -p  /Users/hubq/Downloads/Project/tangEpiPipeline/tangEpiPipelineInstall/out
chmod 777 /Users/hubq/Downloads/Project/tangEpiPipeline/tangEpiPipelineInstall/out

docker run  -v /Users/hubq/Downloads/Project/tangEpiPipeline/tangEpiPipelineInstall/test_fq/:/fastq -v /Users/hubq/Downloads/Project/tangEpiPipeline/tangEpiPipelineInstall/out:/home/analyzer/project -v /Volumes/MacintoshHD/Users/hubq/Downloads/FileZilla/DataBase/mm10/:/home/analyzer/database_ChIP/mm10  -v /Users/hubq/Downloads/Project/tangEpiPipeline/tangEpiPipelineInstall/settings/:/settings/ --env ref=mm10 --env type=ChIP tanginstall:v1

docker run  -v /Users/hubq/Downloads/Project/tangEpiPipeline/tangEpiPipelineInstall/test_fq_RNA/:/fastq -v /Users/hubq/Downloads/Project/tangEpiPipeline/tangEpiPipelineInstall/outRNA:/home/analyzer/project -v /Volumes/MacintoshHD/Users/hubq/Downloads/FileZilla/DataBase/mm10/:/home/analyzer/database_RNA/mm10  -v /Users/hubq/Downloads/Project/tangEpiPipeline/tangEpiPipelineInstall/settings/:/settings/ --env ref=mm10 --env type=RNA tanginstall:v1
```
