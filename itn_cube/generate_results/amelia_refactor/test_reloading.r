

# dsub --provider google-v2 --project my-test-project-210811 --image gcr.io/my-test-project-210811/map_geospatial --regions europe-west1 --label "type=itn_cube" --machine-type n1-standard-1 --logging gs://map_data_z/users/amelia/logs --input-recursive main_dir=gs://map_data_z/users/amelia/itn_cube/predictions --input CODE=gs://map_data_z/users/amelia/itn_cube/code/test_reloading.r --output-recursive output_dir=gs://map_data_z/users/amelia/itn_cube/predictions --command 'Rscript ${CODE}'



library(data.table)

main_dir <- Sys.getenv("main_dir")
output_dir <- Sys.getenv("output_dir")
out_dir <- file.path(output_dir, "test_output_saving.csv")

test_data <- data.table(name=c("Boy", "Snow", "Bird"),
                        chapter=1:3)

write.csv(test_data, out_dir, row.names = F)
print(paste("data saved to", out_dir))
list.files(output_dir)

print("trying to reload data")
new_data <- fread(out_dir)
print("data reloaded successfully")
print(new_data)