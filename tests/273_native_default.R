library(landisR.fire)
old <- options(landisR.fire.backend = NULL)
status <- get(".basefire_native_status_270B", asNamespace("landisR.fire"))
stopifnot(identical(status()$production_fire_backend, "native"),
          isTRUE(status()$fused_native_backend_active), all(status()$compiled_entrypoints))
options(landisR.fire.backend = "legacy")
stopifnot(identical(status()$production_fire_backend, "legacy"),
          !isTRUE(status()$fused_native_backend_active))
options(old)
