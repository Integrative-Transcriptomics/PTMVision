from musialweb import app
from flask import request, session, render_template, send_file
from flask_session import Session
from datetime import timedelta
from dotenv import load_dotenv
import json, zlib, os, subprocess, shutil, numpy, brotli, random, string

""" Load variables from local file system. """
load_dotenv()

""" Set session configuration parameters. """
app.config["SESSION_TYPE"] = "filesystem"
app.config["SESSION_KEY_PREFIX"] = "musial:"
app.config["SESSION_COOKIE_NAME"] = "musial_session"
app.config["SESSION_FILE_DIR"] = "./musialweb/session/"
app.config["SESSION_USE_SIGNER"] = True
app.config["PERMANENT_SESSION_LIFETIME"] = timedelta(days=5.0)
app.config["SESSION_FILE_THRESHOLD"] = 15000
app.config["SECRET_KEY"] = os.getenv("SECRET_KEY")
app.config["MAX_CONTENT_LENGTH"] = 200 * 1024 * 1024  # Limit content lengths to 200 MB.
""" Set API parameters. """
api_parameters = {
    "SUCCESS_CODE": "0",  # Response return code for successful requests.
    "ERROR_CODE": "1",  # Response return code for failed requests (due to server internal error).
    "APPLICATION_ISSUE_CODE": "2",  # Response return code for failed requests (due to application error).
    "RESULT_KEY": os.getenv("RESULT_KEY"),
    "APPLICATION_ERROR_LOG_KEY": os.getenv("APPLICATION_ERROR_LOG_KEY"),
    "SERVER_ERROR_LOG_KEY": os.getenv("SERVER_ERROR_LOG_KEY"),
    "APPLICATION_RUN_LOG_KEY": os.getenv("APPLICATION_RUN_LOG_KEY"),
    "URL": os.getenv("URL"),
}
""" Set DEBUG to true in order to display log as console output. """
DEBUG = bool(int(os.getenv("DEBUG")))
""" Start session """
Session(app)


@app.route("/has_session", methods=["GET"])
def has_session():
    """
    GET route to check whether an active session exists.
    """
    if api_parameters["RESULT_KEY"] in session:
        if session[api_parameters["RESULT_KEY"]] == "":
            return api_parameters["APPLICATION_ISSUE_CODE"]
        else:
            return api_parameters["SUCCESS_CODE"]
    else:
        return api_parameters["ERROR_CODE"]


@app.route("/start_session", methods=["POST"])
def start_session():
    """
    POST route to process user submission.

    The data content of the request has to comply with the MUSIAL run specification (TODO URL), but
    instead of each file path specification the file content is stored. The request's data content is
    expected to be zlib compressed, base64 encoded. The specified file contents will be written to
    the local (server side) file system together with the adjustes MUSIAL run specification.
    MUSIAL is run on the provided data and the generated (TODO variants dicionary) is stored in a opened
    user session, if successfull.
    """
    # Generate unique hex string to use as directory name in the local file system.
    unique_hex_key = generate_random_string()
    # Variables to store output of MUSIAL run.
    stdout = ""
    stderr = ""
    result = ""
    response_code = api_parameters["SUCCESS_CODE"]
    try:
        # Generate directory to store data temporary in the local file system.
        os.mkdir("./musialweb/tmp/" + unique_hex_key)
        # Inflate the request data and transform into python dictionary.
        inflated_request_data = zlib.decompress(request.data)
        json_string_request_data = inflated_request_data.decode("utf8")
        json_request_data = json.loads(json_string_request_data)
        # Write reference .fasta to local file and set path in run specification.
        with open(
            "./musialweb/tmp/" + unique_hex_key + "/reference.fasta", "w+"
        ) as reference_fasta_file:
            reference_fasta_file.write(json_request_data["referenceFASTA"])
            json_request_data["referenceFASTA"] = (
                "./musialweb/tmp/" + unique_hex_key + "/reference.fasta"
            )
        # Write reference .gff3 to local file and set path in run specification.
        with open(
            "./musialweb/tmp/" + unique_hex_key + "/reference.gff3", "w+"
        ) as reference_gff3_file:
            reference_gff3_file.write(json_request_data["referenceGFF"])
            json_request_data["referenceGFF"] = (
                "./musialweb/tmp/" + unique_hex_key + "/reference.gff3"
            )
        # For each specified feature, write .pdb to local file and set path in run specification, if provided.
        for feature in json_request_data["features"].keys():
            if "pdbFile" in json_request_data["features"][feature]:
                with open(
                    "./musialweb/tmp/" + unique_hex_key + "/" + feature + ".pdb", "w+"
                ) as feature_pdb_file:
                    feature_pdb_file.write(
                        json_request_data["features"][feature]["pdbFile"]
                    )
                    json_request_data["features"][feature]["pdbFile"] = (
                        "./musialweb/tmp/" + unique_hex_key + "/" + feature + ".pdb"
                    )
        # For each specified sample, write .vcf to local file and set path in run specification.
        for sample in json_request_data["samples"].keys():
            with open(
                "./musialweb/tmp/" + unique_hex_key + "/" + sample + ".vcf", "w+"
            ) as sample_vcf_file:
                sample_vcf_file.write(json_request_data["samples"][sample]["vcfFile"])
                json_request_data["samples"][sample]["vcfFile"] = (
                    "./musialweb/tmp/" + unique_hex_key + "/" + sample + ".vcf"
                )
        # Write the adjusted request (i.e. used as MUSIAL run specification) to local file.
        with open(
            "./musialweb/tmp/" + unique_hex_key + "/config.json", "w+"
        ) as run_config_file:
            json_request_data["MODULE"] = "BUILD"
            json_request_data["outputFile"] = (
                "./musialweb/tmp/" + unique_hex_key + "/out.vdict.json"
            )
            json.dump([json_request_data], run_config_file)
        # Run MUSIAL on the specified data.
        process = subprocess.Popen(
            [
                os.getenv("JAVA_PATH"),
                "-jar",
                "./musialweb/MUSIAL-v2.1.jar",
                "-c",
                "./musialweb/tmp/" + unique_hex_key + "/config.json",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        stdout, stderr = process.communicate()
        stdout = stdout.decode(encoding="utf-8")
        stderr = stderr.decode(encoding="utf-8")
        # If any error was raised during MUSIAL, set response code to 1 (failed).
        if stderr != "":
            response_code = api_parameters["APPLICATION_ISSUE_CODE"]
            session[api_parameters["APPLICATION_ERROR_LOG_KEY"]] = stderr
            if DEBUG:
                print("\033[41m APPLICATION ERROR \033[0m")
                print(session[api_parameters["APPLICATION_ERROR_LOG_KEY"]])
        # Else, parse and store output of MUSIAL.
        else:
            with open(
                "./musialweb/tmp/" + unique_hex_key + "/out.vdict.json", "r"
            ) as run_out_file:
                json_result_data = json.load(run_out_file)
                result = json_result_data
    # If any error is thrown by the server, set response code to 1 (failed).
    except Exception as e:
        response_code = api_parameters["ERROR_CODE"]
        session[api_parameters["SERVER_ERROR_LOG_KEY"]] = str(e)
        if DEBUG:
            print("\033[41m SERVER ERROR \033[0m")
            print(session[api_parameters["SERVER_ERROR_LOG_KEY"]])
    finally:
        # Remove temporary files, store results and log in session and return response code.
        shutil.rmtree("./musialweb/tmp/" + unique_hex_key)
        session[api_parameters["RESULT_KEY"]] = result
        session[api_parameters["APPLICATION_RUN_LOG_KEY"]] = stdout
        if DEBUG:
            print("\033[46m APPLICATION LOG \033[0m")
            print(session[api_parameters["APPLICATION_RUN_LOG_KEY"]])
        return response_code


@app.route("/log", methods=["GET"])
def get_log():
    log = {}
    if api_parameters["APPLICATION_RUN_LOG_KEY"] in session:
        log[api_parameters["APPLICATION_RUN_LOG_KEY"]] = session[
            api_parameters["APPLICATION_RUN_LOG_KEY"]
        ]
    if api_parameters["APPLICATION_ERROR_LOG_KEY"] in session:
        log[api_parameters["APPLICATION_ERROR_LOG_KEY"]] = session[
            api_parameters["APPLICATION_ERROR_LOG_KEY"]
        ]
    if api_parameters["SERVER_ERROR_LOG_KEY"] in session:
        log[api_parameters["SERVER_ERROR_LOG_KEY"]] = session[
            api_parameters["SERVER_ERROR_LOG_KEY"]
        ]
    if len(log.keys()) > 0:
        return json.dumps(log)
    else:
        return api_parameters["ERROR_CODE"]


@app.route("/get_samples_table", methods=["GET"])
def get_samples_table():
    """
    GET route to return sample information as JSON to display in tabular form.
    """
    if api_parameters["RESULT_KEY"] in session:
        row_data = []
        for sample_identifier in session[api_parameters["RESULT_KEY"]]["samples"]:
            sample_row_data = {
                "Name": session[api_parameters["RESULT_KEY"]]["samples"][
                    sample_identifier
                ]["name"],
                "Substitutions": 0,
                "Deletions": 0,
                "Insertions": 0,
                "Reference Alleles[%]": [0, 0],
                "Reference Proteoforms[%]": [0, 0],
            }
            rec_sub = 0
            rec_del = 0
            rec_ins = 0
            for annotation_field, annotation_value in session[
                api_parameters["RESULT_KEY"]
            ]["samples"][sample_identifier]["annotations"].items():
                if annotation_field.startswith("AL!"):
                    annotation_key = annotation_field.replace("AL!", "AL§")
                    sample_row_data["Reference Alleles[%]"][1] += 1
                    if annotation_value.endswith("REFERENCE"):
                        sample_row_data["Reference Alleles[%]"][0] += 1
                    rec_sub += int(
                        session[api_parameters["RESULT_KEY"]]["features"][
                            annotation_field.split("!")[1]
                        ]["alleles"][annotation_value]["annotations"]["NO_SUB"]
                    )
                    rec_del += int(
                        session[api_parameters["RESULT_KEY"]]["features"][
                            annotation_field.split("!")[1]
                        ]["alleles"][annotation_value]["annotations"]["NO_DEL"]
                    )
                    rec_ins += int(
                        session[api_parameters["RESULT_KEY"]]["features"][
                            annotation_field.split("!")[1]
                        ]["alleles"][annotation_value]["annotations"]["NO_INS"]
                    )
                elif annotation_field.startswith("PF!"):
                    annotation_key = annotation_field.replace("PF!", "PF§")
                    sample_row_data["Reference Proteoforms[%]"][1] += 1
                    if annotation_value.endswith("REFERENCE"):
                        sample_row_data["Reference Proteoforms[%]"][0] += 1
                else:
                    annotation_key = "AN§" + annotation_field
                sample_row_data[annotation_key.replace(".", ":")] = annotation_value
            sample_row_data["Substitutions"] = rec_sub
            sample_row_data["Deletions"] = rec_del
            sample_row_data["Insertions"] = rec_ins
            row_data.append(sample_row_data)
        return json.dumps(row_data)
    else:
        return "1"


@app.route("/get_genes_table", methods=["GET"])
def get_genes_table():
    """
    GET route to return sample information as JSON to display in tabular form.
    """
    if (
        api_parameters["RESULT_KEY"] in session
        and session[api_parameters["RESULT_KEY"]] != ""
    ):
        row_data = []
        for feature_identifier in session[api_parameters["RESULT_KEY"]]["features"]:
            gene_row_data = {
                "Name": session[api_parameters["RESULT_KEY"]]["features"][
                    feature_identifier
                ]["name"],
                "Chromosome": session[api_parameters["RESULT_KEY"]]["features"][
                    feature_identifier
                ]["chromosome"],
                "Start": session[api_parameters["RESULT_KEY"]]["features"][
                    feature_identifier
                ]["start"],
                "End": session[api_parameters["RESULT_KEY"]]["features"][
                    feature_identifier
                ]["end"],
                "Orientation": "+"
                if session[api_parameters["RESULT_KEY"]]["features"][
                    feature_identifier
                ]["isSense"]
                else "-",
                "Coding": "True"
                if session[api_parameters["RESULT_KEY"]]["features"][
                    feature_identifier
                ]["isCodingSequence"]
                else "False",
                "Structure": "True"
                if "structure"
                in session[api_parameters["RESULT_KEY"]]["features"][feature_identifier]
                else "False",
                "Alleles": len(
                    session[api_parameters["RESULT_KEY"]]["features"][
                        feature_identifier
                    ]["alleles"].keys()
                ),
                "Alleles' Variability": [],
                "Proteoforms": len(
                    session[api_parameters["RESULT_KEY"]]["features"][
                        feature_identifier
                    ]["proteoforms"].keys()
                ),
                "Proteoforms' Variability": [],
            }
            per_allele_variability = []
            for allele_object in session[api_parameters["RESULT_KEY"]]["features"][
                feature_identifier
            ]["alleles"].values():
                per_allele_variability.append(
                    float(allele_object["annotations"]["VAR_POS"])
                )
            hist_counts, _ = numpy.histogram(
                per_allele_variability, bins=numpy.linspace(0.0, 10.0, num=21)
            )
            gene_row_data["Alleles' Variability"] = [int(c) for c in list(hist_counts)]
            per_proteoform_variability = []
            for proteoform_object in session[api_parameters["RESULT_KEY"]]["features"][
                feature_identifier
            ]["proteoforms"].values():
                per_proteoform_variability.append(
                    float(proteoform_object["annotations"]["VAR_POS"])
                )
            hist_counts, _ = numpy.histogram(
                per_proteoform_variability, bins=numpy.linspace(0.0, 10.0, num=21)
            )
            gene_row_data["Proteoforms' Variability"] = [
                int(c) for c in list(hist_counts)
            ]
            for annotation_field, annotation_value in session[
                api_parameters["RESULT_KEY"]
            ]["features"][feature_identifier]["annotations"].items():
                gene_row_data[
                    "AN§" + annotation_field.replace(".", ":")
                ] = annotation_value

            row_data.append(gene_row_data)
        return json.dumps(row_data)
    else:
        return api_parameters["ERROR_CODE"]


@app.route("/get_variants_table", methods=["GET"])
def get_variants_table():
    """
    GET route to return variants information as JSON to display in tabular form.
    """
    if (
        api_parameters["RESULT_KEY"] in session
        and session[api_parameters["RESULT_KEY"]] != ""
    ):
        row_data = []
        for position in session[api_parameters["RESULT_KEY"]]["nucleotideVariants"]:
            for content in session[api_parameters["RESULT_KEY"]]["nucleotideVariants"][
                position
            ]:
                alt = str(content)
                ref = str(
                    session[api_parameters["RESULT_KEY"]]["nucleotideVariants"][
                        position
                    ][content]["annotations"]["REF_CONT"]
                )
                type = "N/A"
                if len(alt) == 1 == len(ref):
                    type = "Substitution"
                elif len(alt.replace("-", "")) < len(ref):
                    type = "Deletion"
                else:
                    type = "Insertion"
                variant_row_data = {
                    "Position": position,
                    "Type": type,
                    "Reference": ref,
                    "Alternative": alt,
                    "Frequency[%]": session[api_parameters["RESULT_KEY"]][
                        "nucleotideVariants"
                    ][position][content]["annotations"]["FRQ"],
                }
                effects = set([])
                occurences = set([])
                for affected_feature in session[api_parameters["RESULT_KEY"]][
                    "nucleotideVariants"
                ][position][content]["occurrence"]:
                    feature, allele = affected_feature.split("!")
                    effects.add(feature)
                    occurences.add(allele)
                variant_row_data["Effect"] = ";".join(effects)
                variant_row_data["Occurrence"] = ";".join(occurences)
                for annotation_field, annotation_value in session[
                    api_parameters["RESULT_KEY"]
                ]["nucleotideVariants"][position][content]["annotations"].items():
                    if annotation_field.startswith("SNPEFF_"):
                        variant_row_data[
                            annotation_field.replace(".", ":")
                        ] = annotation_value
                row_data.append(variant_row_data)
        return json.dumps(row_data)
    else:
        return api_parameters["ERROR_CODE"]


@app.route("/get_feature_annotation", methods=["GET"])
def get_feature_annotation():
    """
    GET route to return annotation of a features proteoform or allele.
    """
    if (
        api_parameters["RESULT_KEY"] in session
        and session[api_parameters["RESULT_KEY"]] != ""
    ):
        data = {}
        feature, target = request.args.get("target").split("%")
        feature = feature.replace(":", ".")
        if target.startswith("AL"):
            target_data = session[api_parameters["RESULT_KEY"]]["features"][feature][
                "alleles"
            ][target]
            data = {
                "Name": target_data["name"],
                "Deletions": target_data["annotations"]["NO_DEL"],
                "Substitutions": target_data["annotations"]["NO_SUB"],
                "Insertions": target_data["annotations"]["NO_INS"],
                "Frequency[%]": target_data["annotations"]["FRQ"],
                "Variability[%]": target_data["annotations"]["VAR_POS"],
                "Variants": target_data["annotations"]["VARIANTS"],
            }
        elif target.startswith("PF"):
            target_data = session[api_parameters["RESULT_KEY"]]["features"][feature][
                "proteoforms"
            ][target]
            data = {
                "Name": target_data["name"],
                "Deletions": target_data["annotations"]["NO_DEL"],
                "Substitutions": target_data["annotations"]["NO_SUB"],
                "Insertions": target_data["annotations"]["NO_INS"],
                "Frequency [%]": target_data["annotations"]["FRQ"],
                "Variability [%]": target_data["annotations"]["VAR_POS"],
                "Conglomeration [p-value]": target_data["annotations"]["CONG_IDX"],
                "Novel Termination": target_data["annotations"]["DIV_TERM_POS"],
                "Truncation [%]": target_data["annotations"]["DIV_TERM_TRC_PRC"],
                "Variants": target_data["annotations"]["VARIANTS"],
            }
        return json.dumps(data)
    else:
        return api_parameters["ERROR_CODE"]


@app.route("/protein_dashboard", methods=["GET"])
def get_protein_dashboard():
    target_identifier = request.args.get("target")
    if (
        api_parameters["RESULT_KEY"] in session
        and session[api_parameters["RESULT_KEY"]] != ""
        and target_identifier in session[api_parameters["RESULT_KEY"]]["features"]
    ):
        target_data = session[api_parameters["RESULT_KEY"]]["features"][
            target_identifier
        ]
        return render_template(
            "protein_dashboard.html", data=target_data, target=target_identifier
        )
    else:
        return api_parameters["ERROR_CODE"]


@app.route("/result", methods=["GET"])
def get_result():
    # Generate unique hex string to use as directory name in the local file system.
    unique_hex_key = generate_random_string()
    response_code = api_parameters["SUCCESS_CODE"]
    try:
        # Generate directory to store data temporary in the local file system.
        os.mkdir("./musialweb/tmp/" + unique_hex_key)
        with open(
            "./musialweb/tmp/" + unique_hex_key + "/result.json.br", "wb+"
        ) as compressed_result:
            # Write compressed session result.
            compressed_result.write(
                brotli.compress(
                    json.dumps(session[api_parameters["RESULT_KEY"]]).encode("utf-8")
                )
            )
            return send_file(
                "./tmp/" + unique_hex_key + "/result.json.br",
                as_attachment=True,
            )
    except Exception as e:
        response_code = api_parameters["ERROR_CODE"]
        session[api_parameters["SERVER_ERROR_LOG_KEY"]] = str(e)
        if DEBUG:
            print("\033[41m SERVER ERROR \033[0m")
            print(session[api_parameters["SERVER_ERROR_LOG_KEY"]])
        return response_code
    finally:
        shutil.rmtree("./musialweb/tmp/" + unique_hex_key)


@app.route("/download_sequences", methods=["POST"])
def get_sequences():
    # Generate unique hex string to use as directory name in the local file system.
    unique_hex_key = generate_random_string()
    # Variables to store output of MUSIAL run.
    stdout = ""
    stderr = ""
    response_code = api_parameters["SUCCESS_CODE"]
    try:
        # Generate directory to store data temporary in the local file system.
        os.mkdir("./musialweb/tmp/" + unique_hex_key)
        # Inflate the request data and transform into python dictionary.
        inflated_request_data = zlib.decompress(request.data)
        json_string_request_data = inflated_request_data.decode("utf8")
        json_request_data = json.loads(json_string_request_data)
        # Write session result to temporary file in the local file system.
        with open(
            "./musialweb/tmp/" + unique_hex_key + "/result.vdict.json", "w+"
        ) as tmp_result_file:
            json.dump(session[api_parameters["RESULT_KEY"]], tmp_result_file)
            json_request_data["inputFile"] = (
                "./musialweb/tmp/" + unique_hex_key + "/result.vdict.json"
            )
        # Set output file to temporary file path in local file system.
        os.mkdir("./musialweb/tmp/" + unique_hex_key + "/sequences")
        json_request_data["outputDirectory"] = (
            "./musialweb/tmp/" + unique_hex_key + "/sequences/"
        )
        if (
            "includeSamples" in json_request_data
            and len(json_request_data["includeSamples"]) > 0
        ):
            json_request_data["samples"] = json_request_data["includeSamples"]
            json_request_data.pop("includeSamples")
        elif (
            "excludeSamples" in json_request_data
            and len(json_request_data["excludeSamples"]) > 0
        ):
            samples_to_include = list(
                session[api_parameters["RESULT_KEY"]]["samples"].keys()
            )
            for sample_to_exclude in json_request_data["excludeSamples"]:
                samples_to_include.remove(sample_to_exclude)
            json_request_data["samples"] = samples_to_include
            json_request_data.pop("excludeSamples")
        if (
            "includeFeatures" in json_request_data
            and len(json_request_data["includeFeatures"]) > 0
        ):
            json_request_data["features"] = json_request_data["includeFeatures"]
            json_request_data.pop("includeFeatures")
        elif (
            "excludeFeatures" in json_request_data
            and len(json_request_data["excludeFeatures"]) > 0
        ):
            features_to_include = list(
                session[api_parameters["RESULT_KEY"]]["features"].keys()
            )
            for feature_to_exclude in json_request_data["excludeFeatures"]:
                features_to_include.remove(feature_to_exclude)
            json_request_data["features"] = features_to_include
            json_request_data.pop("excludeFeatures")
        # Write the adjusted request (i.e. used as MUSIAL run specification) to local file.
        with open(
            "./musialweb/tmp/" + unique_hex_key + "/config.json", "w+"
        ) as run_config_file:
            json.dump([json_request_data], run_config_file)
        # Run MUSIAL on the specified data.
        process = subprocess.Popen(
            [
                os.getenv("JAVA_PATH"),
                "-jar",
                "./musialweb/MUSIAL-v2.1.jar",
                "-c",
                "./musialweb/tmp/" + unique_hex_key + "/config.json",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        stdout, stderr = process.communicate()
        stdout = stdout.decode(encoding="utf-8")
        stderr = stderr.decode(encoding="utf-8")
        # If any error was raised during MUSIAL, set response code to 1 (failed).
        if stderr != "":
            response_code = api_parameters["APPLICATION_ISSUE_CODE"]
            session[api_parameters["APPLICATION_ERROR_LOG_KEY"]] = stderr
            if DEBUG:
                print("\033[41m APPLICATION ERROR \033[0m")
                print(session[api_parameters["APPLICATION_ERROR_LOG_KEY"]])
            return response_code
        # Else, compress output and send to client.
        else:
            shutil.make_archive(
                "./musialweb/tmp/" + unique_hex_key + "/sequences",
                "zip",
                "./musialweb/tmp/" + unique_hex_key + "/sequences/",
            )
            return send_file(
                "./tmp/" + unique_hex_key + "/sequences.zip",
                as_attachment=True,
                mimetype="application/octet-stream",
            )
    # If any error is thrown by the server, set response code to 1 (failed).
    except Exception as e:
        response_code = api_parameters["ERROR_CODE"]
        session[api_parameters["SERVER_ERROR_LOG_KEY"]] = str(e)
        if DEBUG:
            print("\033[41m SERVER ERROR \033[0m")
            print(session[api_parameters["SERVER_ERROR_LOG_KEY"]])
        return response_code
    finally:
        # Remove temporary files, store results and log in session and return response code.
        shutil.rmtree("./musialweb/tmp/" + unique_hex_key)
        session[api_parameters["APPLICATION_RUN_LOG_KEY"]] = stdout
        if DEBUG:
            print("\033[46m APPLICATION LOG \033[0m")
            print(session[api_parameters["APPLICATION_RUN_LOG_KEY"]])


@app.route("/download_tables", methods=["POST"])
def get_tables():
    # Generate unique hex string to use as directory name in the local file system.
    unique_hex_key = generate_random_string()
    # Variables to store output of MUSIAL run.
    stdout = ""
    stderr = ""
    response_code = api_parameters["SUCCESS_CODE"]
    try:
        # Generate directory to store data temporary in the local file system.
        os.mkdir("./musialweb/tmp/" + unique_hex_key)
        # Inflate the request data and transform into python dictionary.
        inflated_request_data = zlib.decompress(request.data)
        json_string_request_data = inflated_request_data.decode("utf8")
        json_request_data = json.loads(json_string_request_data)
        # Write session result to temporary file in the local file system.
        with open(
            "./musialweb/tmp/" + unique_hex_key + "/result.vdict.json", "w+"
        ) as tmp_result_file:
            json.dump(session[api_parameters["RESULT_KEY"]], tmp_result_file)
            json_request_data["inputFile"] = (
                "./musialweb/tmp/" + unique_hex_key + "/result.vdict.json"
            )
        # Set output file to temporary file path in local file system.
        os.mkdir("./musialweb/tmp/" + unique_hex_key + "/tables")
        json_request_data["outputDirectory"] = (
            "./musialweb/tmp/" + unique_hex_key + "/tables/"
        )
        if (
            "includeSamples" in json_request_data
            and len(json_request_data["includeSamples"]) > 0
        ):
            json_request_data["samples"] = json_request_data["includeSamples"]
            json_request_data.pop("includeSamples")
        elif (
            "excludeSamples" in json_request_data
            and len(json_request_data["excludeSamples"]) > 0
        ):
            samples_to_include = list(
                session[api_parameters["RESULT_KEY"]]["samples"].keys()
            )
            for sample_to_exclude in json_request_data["excludeSamples"]:
                samples_to_include.remove(sample_to_exclude)
            json_request_data["samples"] = samples_to_include
            json_request_data.pop("excludeSamples")
        if (
            "includeFeatures" in json_request_data
            and len(json_request_data["includeFeatures"]) > 0
        ):
            json_request_data["features"] = json_request_data["includeFeatures"]
            json_request_data.pop("includeFeatures")
        elif (
            "excludeFeatures" in json_request_data
            and len(json_request_data["excludeFeatures"]) > 0
        ):
            features_to_include = list(
                session[api_parameters["RESULT_KEY"]]["features"].keys()
            )
            for feature_to_exclude in json_request_data["excludeFeatures"]:
                features_to_include.remove(feature_to_exclude)
            json_request_data["features"] = features_to_include
            json_request_data.pop("excludeFeatures")
        # Write the adjusted request (i.e. used as MUSIAL run specification) to local file.
        with open(
            "./musialweb/tmp/" + unique_hex_key + "/config.json", "w+"
        ) as run_config_file:
            json.dump([json_request_data], run_config_file)
        # Run MUSIAL on the specified data.
        process = subprocess.Popen(
            [
                os.getenv("JAVA_PATH"),
                "-jar",
                "./musialweb/MUSIAL-v2.1.jar",
                "-c",
                "./musialweb/tmp/" + unique_hex_key + "/config.json",
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        stdout, stderr = process.communicate()
        stdout = stdout.decode(encoding="utf-8")
        stderr = stderr.decode(encoding="utf-8")
        # If any error was raised during MUSIAL, set response code to 1 (failed).
        if stderr != "":
            response_code = api_parameters["APPLICATION_ISSUE_CODE"]
            session[api_parameters["APPLICATION_ERROR_LOG_KEY"]] = stderr
            if DEBUG:
                print("\033[41m APPLICATION ERROR \033[0m")
                print(session[api_parameters["APPLICATION_ERROR_LOG_KEY"]])
            return response_code
        # Else, compress output and send to client.
        else:
            shutil.make_archive(
                "./musialweb/tmp/" + unique_hex_key + "/tables",
                "zip",
                "./musialweb/tmp/" + unique_hex_key + "/tables/",
            )
            return send_file(
                "./tmp/" + unique_hex_key + "/tables.zip",
                as_attachment=True,
                mimetype="application/octet-stream",
            )
    # If any error is thrown by the server, set response code to 1 (failed).
    except Exception as e:
        response_code = api_parameters["ERROR_CODE"]
        session[api_parameters["SERVER_ERROR_LOG_KEY"]] = str(e)
        if DEBUG:
            print("\033[41m SERVER ERROR \033[0m")
            print(session[api_parameters["SERVER_ERROR_LOG_KEY"]])
        return response_code
    finally:
        # Remove temporary files, store results and log in session and return response code.
        shutil.rmtree("./musialweb/tmp/" + unique_hex_key)
        session[api_parameters["APPLICATION_RUN_LOG_KEY"]] = stdout
        if DEBUG:
            print("\033[46m APPLICATION LOG \033[0m")
            print(session[api_parameters["APPLICATION_RUN_LOG_KEY"]])


@app.route("/example_data", methods=["GET"])
def get_example_data():
    return send_file("./static/resources/example_data.zip", as_attachment=True)


@app.route("/example_session", methods=["GET"])
def get_example_session():
    """
    GET route to start an example session.
    """
    # Variables to store output of MUSIAL run.
    result = ""
    response_code = api_parameters["SUCCESS_CODE"]
    try:
        # Load static example session.
        with open(
            "./musialweb/static/resources/example_session.json", "r"
        ) as example_session_file:
            result = json.load(example_session_file)
    # If any error is thrown by the server, set response code to 1 (failed).
    except Exception as e:
        response_code = api_parameters["ERROR_CODE"]
        session[api_parameters["SERVER_ERROR_LOG_KEY"]] = str(e)
        if DEBUG:
            print("\033[41m SERVER ERROR \033[0m")
            print(session[api_parameters["SERVER_ERROR_LOG_KEY"]])
    finally:
        session[api_parameters["RESULT_KEY"]] = result
        return response_code


def generate_random_string():
    return "".join(
        random.SystemRandom().choice(string.ascii_letters + string.digits)
        for _ in range(10)
    )
