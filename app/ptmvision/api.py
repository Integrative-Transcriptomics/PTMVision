from ptmvision import app
from ptmvision import utils
from ptmvision.uniprot_id_mapping import (
    submit_id_mapping,
    check_id_mapping_results_ready,
    get_id_mapping_results_link,
    get_id_mapping_results_search,
)
from flask import request, session
from flask_session import Session
from cachelib import FileSystemCache
from dotenv import load_dotenv
from io import StringIO
from itertools import combinations
from copy import deepcopy
from datetime import datetime
import json, zlib, os, base64, traceback

""" Load variables from local file system. """
load_dotenv()

""" Definition of session keys """
MODIFICATIONS_DATA = "bW9kaWZpY2F0aW9uc19kYXRhX2ZyYW1l"
BASEPATH = "./app/ptmvision" # Use this for local development.
#BASEPATH = "/app/ptmvision" # Use this for deployment.

""" Set session configuration parameters. """
app.config["SESSION_COOKIE_NAME"] = "PTMVision"
app.config["SESSION_TYPE"] = "cachelib"
app.config["SESSION_PERMANENT"] = False
app.config["SESSION_USE_SIGNER"] = True
app.config["SECRET_KEY"] = os.getenv("SIGNER")
app.config["SESSION_CACHELIB"] = FileSystemCache(cache_dir = BASEPATH + '/session/', threshold=1000)
app.config["SESSION_CLEANUP_N_REQUESTS"] = 100
app.config["MAX_CONTENT_LENGTH"] = 100 * 1024 * 1024  # Limit content lengths to 100 MB.

""" Set API parameters. """
DEBUG = os.getenv("DEBUG")

""" Start session """
Session(app)

@app.route("/example_session", methods=["GET"])
def example_session():
    try :
        session.clear( )
        with open( BASEPATH + "/static/resources/example_modifications_data.json", "r" ) as example_modifications_data :
            session[MODIFICATIONS_DATA] = json.load( example_modifications_data )
        return "Ok", 200
    except Exception as e :
        return "Failed request '/example_session': " + _format_exception(e), 500


@app.route("/download_session", methods=["GET"])
def download_session():
    try :
        zlib_compress = zlib.compressobj( 6, zlib.DEFLATED, zlib.MAX_WBITS )
        compressed_session_bytes = zlib_compress.compress( bytes(json.dumps(session[MODIFICATIONS_DATA]), "utf-8") ) + zlib_compress.flush( )
        encoded_session = base64.b64encode( compressed_session_bytes ).decode("ascii")
        return encoded_session, 200
    except Exception as e :
        return "Failed request '/download_session': " + _format_exception(e), 500


@app.route("/restart_session", methods=["POST"])
def restart_session():
    try :
        session.clear( )
        session_data = zlib.decompress( base64.b64decode( request.data ) ).decode( )
        session[MODIFICATIONS_DATA] = json.loads(session_data)
        return "Ok", 200
    except Exception as e :
        return "Failed request '/example_session': " + _format_exception(e), 500


@app.route("/process_search_engine_output", methods=["POST"])
def process_search_engine_output():
    """Route to process CSV data or SAGE, ionbot or MSFragger search engine output."""
    try:
        session.clear( )
        json_request_data = _request_to_json(request.data)
        json_user_data = utils.parse_user_input(
            StringIO(json_request_data["content"]),
            json_request_data["contentType"],
            json_request_data["massShiftTolerance"],
            json_request_data["excludeClasses"]
        )
        session[MODIFICATIONS_DATA] = json_user_data
        if DEBUG :
            with open( "./dump.json", "w+" ) as dumpfile :
                dumpfile.write( json.dumps( session[MODIFICATIONS_DATA], indent = 3 ) )
        return "Ok", 200
    except Exception as e:
        return "Failed request '/process_search_engine_output': " + _format_exception(e), 500


@app.route("/get_available_proteins", methods=["GET"])
def get_available_proteins():
    """Route to retrieve all available proteins of the session as a JSON."""
    try :
        protein_entries = []
        if MODIFICATIONS_DATA in session:
            protein_identifiers_unannotated = list(
                filter(
                    lambda _ : "annotation" not in session[MODIFICATIONS_DATA]["proteins"][ _ ],
                    list( session[MODIFICATIONS_DATA]["proteins"].keys( ) )
                )
            )
            if len( protein_identifiers_unannotated ) > 0 :
                annotations = _map_uniprot_identifiers(
                    protein_identifiers_unannotated,
                    "UniProtKB"
                )

                # remove IDs that were demerged in UniProt into multiple
                # we don't know the original sequence, site mapping would go wrong
                seen = set()
                duplicates = set(item["from"] for item in annotations["results"] if item["from"] in seen or seen.add(item["from"]))

                # remove them from results 
                annotations["results"] = [item for item in annotations["results"] if item["from"] not in duplicates]

                if "failedIds" in annotations :
                    annotations["failedIds"].extend(list(duplicates))
                else :
                    annotations["failedIds"] = list(duplicates)
                
                # Number of demerged proteins / failed IDs to show user
                #if env_parameters["DEBUG"] :
                #n_demerged = len(duplicates)
                #n_failed = len(annotations["failedIds"])

                if "results" in annotations :
                    for entry in annotations[ "results" ] :
                        _add_annotation_to_protein( entry["from"], entry["to"] )
                if "failedIds" in annotations :
                    for failed_id in annotations[ "failedIds" ] :
                        _add_annotation_to_protein( failed_id, { }, failed = True )

            # Construct entry per protein in input data.
            for protein_identifier in list( session[MODIFICATIONS_DATA]["proteins"].keys( ) ):
                protein_name = _get_protein_name(session[MODIFICATIONS_DATA]["proteins"][protein_identifier]["annotation"])
                protein_length = _get_protein_length(session[MODIFICATIONS_DATA]["proteins"][protein_identifier]["annotation"])
                position_modification_data = session[MODIFICATIONS_DATA]["proteins"][
                    protein_identifier
                ]["positions"]
                modified_positions = []
                modifications = []
                for modified_position in position_modification_data:
                    modified_positions.append(int(modified_position))
                    for modification in position_modification_data[modified_position][
                        "modifications"
                    ]:
                        modifications.append(modification["display_name"])
                modified_positions = sorted(modified_positions)
                modifications = list(set(modifications))
                protein_entry = {
                    "id": protein_identifier,
                    "name": protein_name,
                    "length": protein_length,
                    "modified_positions": len(modified_positions),
                    "unique_modifications": len(modifications),
                    "modifications": "$".join(modifications),
                }
                protein_entries.append(protein_entry)
            if DEBUG :
                    with open( "./dump.json", "w+" ) as dumpfile :
                        dumpfile.write( json.dumps( session[MODIFICATIONS_DATA], indent = 3 ) )
            return protein_entries, 200
        else :
            raise Exception("Faulty session data.")
    except Exception as e:
        return "Failed request '/get_available_proteins': " + _format_exception(e), 500


@app.route("/get_overview_data", methods=["GET"])
def get_overview_data():
    try :
        modifications =  { } # Stores all present modifications together with their count.
        modification_co_occurrence = { } # Maps pairs of modification display names to their co-occurrence count.
        modification_classification_counts = { }
        if MODIFICATIONS_DATA in session:
            protein_identifiers = [_ for _ in session[MODIFICATIONS_DATA]["proteins"]]
            dataset_size = len( protein_identifiers )
            for protein_identifier in protein_identifiers:
                modification_data = session[MODIFICATIONS_DATA]["proteins"][
                    protein_identifier
                ]["positions"]
                for position in modification_data:
                    modifications_at_position = []
                    for modification in modification_data[position][
                        "modifications"
                    ]:
                        modification_name = modification["display_name"]
                        modifications.setdefault( modification_name, deepcopy( modification ) )

                        modifications[ modification_name ].setdefault( "count", 0 )
                        modifications[ modification_name ][ "count" ] += 1

                        modifications[ modification_name ].setdefault( "occurrence", [ ] )
                        modifications[ modification_name ][ "occurrence" ].append( protein_identifier )

                        modifications_at_position.append(modification_name)

                        modification_classification = modification[ "modification_classification" ]
                        modification_classification_counts.setdefault( modification_classification, 0 )
                        modification_classification_counts[ modification_classification ] += 1

                    for modifications_pair_tuple in combinations(set(modifications_at_position), 2):
                        modifications_pair = "@".join( sorted( list( modifications_pair_tuple ) ) )
                        modification_co_occurrence.setdefault( modifications_pair, 0 )
                        modification_co_occurrence[ modifications_pair ] += 1

            for modification_name, modification in modifications.items( ) :
                modification[ "frequency" ] = round( ( len( set( modifications[ modification_name ][ "occurrence" ] ) ) / dataset_size ) * 100, 2 )
            del modifications[ modification_name ][ "occurrence" ]

            modifications_by_mass_shift = sorted( list( modifications.keys( ) ), key = lambda k : -modifications[k][ "mass_shift" ] if type( modifications[k][ "mass_shift" ] ) != str else 0.0 )
            modifications_by_count = sorted( list( modifications.keys( ) ), key = lambda k : -modifications[k][ "count" ] )

            return [
                modifications,
                [ modifications_by_mass_shift, modifications_by_count ],
                modification_co_occurrence,
                modification_classification_counts,
                session[MODIFICATIONS_DATA][ "meta_data" ]
            ], 200
        else :
            raise Exception("Faulty session data.")
    except Exception as e:
        return "Failed request '/get_overview_data': " + _format_exception(e), 500


@app.route("/get_protein_data", methods=["POST"])
def get_protein_data():
    try :
        json_request_data = _request_to_json(request.data)
        # Try to fetch PDB format structure for UniProt identifier from AlphaFold database.
        if json_request_data["pdb_text"] != None:
            structure, pdb_text = utils.parse_structure(json_request_data["pdb_text"])
        else:
            if "structure" in session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_pa"] ] :
                pdb_text = utils._brotli_decompress( session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_pa"] ]["structure"] )
                structure = utils.parse_structure( pdb_text )
            else :
                structure, pdb_text = utils.get_structure(json_request_data["uniprot_pa"])
                # Extract protein sequence from structure and store it in session data.
                session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_pa"] ]["sequence"] = utils.get_sequence_from_structure(structure)
                session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_pa"] ][ "structure" ] = utils._brotly_compress(pdb_text)
        if structure != None:
            # Extract annotations for protein from UniProt.
            # NOTE: Temp. deprecated code segment; Data is (experimental) queried from UniProt by initial analysis for all proteins. This may be re-used if performance issues occur.
            # if not "annotation" in session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_pa"] ] :
            #    annotation = _map_uniprot_identifiers(
            #        [json_request_data["uniprot_pa"]], "UniProtKB"
            #    )
            #    session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_pa"] ]["annotation"] = { }
            #    session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_pa"] ]["annotation"] = annotation[ "results" ][ 0 ][ "to" ]
            # Compute contacts from structure and store them in session data.
            if not "contacts" in session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_pa"] ] :
                session[MODIFICATIONS_DATA]["meta_data"]["distance_cutoff"] = int(json_request_data["cutoff"])
                session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_pa"] ][ "contacts" ] = { }
                distance_matrix = utils.get_distance_matrix(structure)
                contacts = utils.get_contacts(
                    distance_matrix,
                    float(json_request_data["cutoff"]),
                )
                for source_position, contacts_list in contacts.items( ) :
                    session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_pa"] ][ "contacts" ][ source_position ] = [ ]
                    for contact_position in contacts_list :
                        session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_pa"] ][ "contacts" ][ source_position ].append( ( contact_position, distance_matrix[ source_position, contact_position ] ) )
            # Construct response.
            response = deepcopy( session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_pa"] ] )
            response[ "structure" ] = utils._brotli_decompress( response[ "structure" ] )
            response[ "meta" ] = session[MODIFICATIONS_DATA]["meta_data"]
            return response, 200
        else :
            return "Error in request '/get_protein_data': No protein structure available.", 303
    except Exception as e:
        return "Failed request '/get_protein_data': " + _format_exception(e), 500


def _request_to_json(request_content: bytes) -> object:
    """
    Decompress zlib encoded JSON request into Python object.

    Uses UTF-8 decoding.

    Parameters:
        request_content (bytes): Zlib encoded JSON request content.

    Returns:
        object: Python object representing the content of the JSON request.
    """
    inflated_request_data = zlib.decompress(request_content)
    json_string_request_data = inflated_request_data.decode("utf8")
    return json.loads(json_string_request_data)


def _map_uniprot_identifiers(identifiers: list, target_db: str) -> dict:
    """
    Map UniProt identifiers to a target database.

    Parameters:
        identifiers (list): List of UniProt identifiers.
        target_db (str): Target database for mapping.

    Returns:
        dict: Results of the identifier mapping.
    """
    if DEBUG :
        print( "\33[33mSTATUS [" + str(datetime.now()) + "]\33[0m\tRetrieve UniProt information ... ", end = '', flush = True)
    job_id = submit_id_mapping(
        from_db="UniProtKB_AC-ID", to_db=target_db, ids=identifiers
    )
    if check_id_mapping_results_ready(job_id):
        link = get_id_mapping_results_link(job_id)
        results = get_id_mapping_results_search(link)
    if DEBUG :
        print( "DONE", end = '', flush = True)


    return results


def _rename_dictionary_entry(D: dict, key_old: str, key_new: str):
    """
    Rename a key in a dictionary.

    Parameters:
        D (dict): The dictionary containing the entry to be renamed.
        key_old (str): The old key to be replaced.
        key_new (str): The new key to replace the old key.
    """
    if key_old != key_new:
        D[key_new] = D[key_old]
        del D[key_old]


def _add_annotation_to_protein(protein_identifier: str, annotation: dict, failed: bool = False):
    """
    Add annotation to a protein entry.

    Parameters:
        protein_identifier (str): Identifier of the protein.
        annotation (dict): Annotation information for the protein.
        failed (bool, optional): Flag indicating if annotation failed. Defaults to False.
    """
    if failed:
        session[MODIFICATIONS_DATA]["proteins"][protein_identifier]["annotation"] = None
    else:
        if protein_identifier != annotation["primaryAccession"]:
            _rename_dictionary_entry(session[MODIFICATIONS_DATA]["proteins"], protein_identifier, annotation["primaryAccession"])
            protein_identifier = annotation["primaryAccession"]
        # Adjust annotation.
        annotation.pop("entryType", None)
        annotation.pop("secondaryAccessions", None)
        annotation.pop("keywords", None)
        annotation.pop("references", None)
        annotation.pop("uniProtKBCrossReferences", None)
        annotation.pop("extraAttributes", None)
        # Set annotation to protein entry.
        session[MODIFICATIONS_DATA]["proteins"][protein_identifier]["annotation"] = annotation


def _get_protein_name(annotation: dict) -> str:
    """
    Get the name of a protein from its annotation.

    Parameters:
        annotation (dict): Annotation information for the protein.

    Returns:
        str: The name of the protein.
    """
    na = "N/A"
    if annotation == None:
        return na
    else:
        if "proteinDescription" in annotation:
            if "recommendedName" in annotation["proteinDescription"]:
                return annotation["proteinDescription"]["recommendedName"]["fullName"]["value"]
            elif "alternativeNames" in annotation["proteinDescription"]:
                return annotation["proteinDescription"]["alternativeNames"][0]["fullName"]["value"]
            elif "submissionNames" in annotation["proteinDescription"]:
                return annotation["proteinDescription"]["submissionNames"][0]["fullName"]["value"]
            else:
                return na
        else:
            return na
        

def _get_protein_length(annotation: dict) -> int:
    """
    Get the length of a protein from its annotation.

    Parameters:
        annotation (dict): Annotation information for the protein.

    Returns:
        int: The length of the protein.
    """
    na = "N/A"
    if annotation == None:
        return na
    else :
        if "sequence" in annotation:
            if "length" in annotation["sequence"] :
                return annotation["sequence"]["length"]
            else :
                return na
        else :
            return na


def _format_exception(e: str) -> str:
    """
    Formats the given exception into a string.

    Args:
        e (str): The exception to be formatted.

    Returns:
        str: The formatted exception as a string.

    Note:
        This function utilizes the global variable `api_parameters` to determine whether to print the exception traceback.
        If `api_parameters["DEBUG"]` is True, the traceback will be printed with color highlighting; otherwise, only the exception message will be formatted.
    """
    if DEBUG:
        # Print the exception traceback with color highlighting
        print("\u001b[31m" + "".join(traceback.format_exception(e)) + "\u001b[0m")
    
    # Format only the exception message
    return "".join(traceback.format_exception_only(e)).strip()