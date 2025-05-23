import os

from conan import ConanFile
from conan.tools.cmake import CMakeToolchain, CMake, cmake_layout
from conan.tools.build import can_run
from conan.tools.files import save
from conan.tools.scm import Git


class HyhoundRecipe(ConanFile):
    name = "hyhound"
    version = "1.0.0"

    # Optional metadata
    license = "LGPLv3"
    author = "Pieter P <pieter.p.dev@outlook.com>"
    url = "https://github.com/kul-optec/hyhound"
    description = "Hyperbolic Householder transformations for Cholesky factorization up- and downdates."
    topics = "scientific software"

    # Binary configuration
    package_type = "library"
    settings = "os", "compiler", "build_type", "arch"
    bool_hyhound_options = {
        "with_ocp": False,
        "with_benchmarks": False,
    }
    options = {
        "shared": [True, False],
        "fPIC": [True, False],
        "real_type": ["double;float", "float;double", "double", "float"],
    } | {k: [True, False] for k in bool_hyhound_options}
    default_options = {
        "shared": False,
        "fPIC": True,
        "real_type": "double;float",
    } | bool_hyhound_options

    # Sources are located in the same place as this recipe, copy them to the recipe
    exports_sources = (
        "CMakeLists.txt",
        "src/*",
        "cmake/*",
        "test/*",
        "benchmarks/*",
        "LICENSE",
        "README.md",
    )

    def export_sources(self):
        git = Git(self)
        status_cmd = "status . --short --no-branch --untracked-files=no"
        dirty = bool(git.run(status_cmd).strip())
        hash = git.get_commit() + ("-dirty" if dirty else "")
        print("Commit hash:", hash)
        save(self, os.path.join(self.export_sources_folder, "commit.txt"), hash)

    generators = ("CMakeDeps",)

    def requirements(self):
        self.requires(
            "guanaqo/1.0.0-alpha.12", transitive_headers=True, transitive_libs=True
        )
        if self.options.with_ocp:
            self.requires("eigen/3.4.0", transitive_headers=True)
        elif self.options.with_benchmarks:
            self.requires("eigen/3.4.0")
        else:
            self.test_requires("eigen/3.4.0")
        if self.options.with_benchmarks:
            self.requires("benchmark/1.8.4")

    def build_requirements(self):
        self.test_requires("gtest/1.15.0")
        self.tool_requires("cmake/[>=3.24 <4]")

    def config_options(self):
        if self.settings.get_safe("os") == "Windows":
            self.options.rm_safe("fPIC")

    def configure(self):
        self.options["guanaqo/*"].with_blas = True

    def layout(self):
        cmake_layout(self)
        self.cpp.build.builddirs.append("")

    def generate(self):
        tc = CMakeToolchain(self)
        for k in self.bool_hyhound_options:
            value = self.options.get_safe(k)
            if value is not None and value.value is not None:
                tc.variables["HYHOUND_" + k.upper()] = bool(value)
        guanaqo = self.dependencies["guanaqo"]
        index_type = guanaqo.options.get_safe("blas_index_type", default="int")
        real_type = str(self.options.real_type)
        print("index_type:", index_type)
        print("real_type: ", real_type)
        tc.variables["HYHOUND_DENSE_INDEX_TYPE"] = index_type
        tc.variables["HYHOUND_DENSE_REAL_TYPE"] = real_type
        if can_run(self):
            tc.variables["HYHOUND_FORCE_TEST_DISCOVERY"] = True
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()
        cmake.test()

    def package(self):
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        self.cpp_info.set_property("cmake_find_mode", "none")
        self.cpp_info.builddirs.append(os.path.join("lib", "cmake", "hyhound"))
