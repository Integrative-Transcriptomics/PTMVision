from ptmvis import app
from ptmvis.backend import utils
from ptmvis.backend.uniprot_id_mapping import (
    submit_id_mapping,
    check_id_mapping_results_ready,
    get_id_mapping_results_link,
    get_id_mapping_results_search,
)
from flask import request, session, send_file
from flask_session import Session
from dotenv import load_dotenv
from io import StringIO
from itertools import combinations
from math import ceil
from copy import deepcopy
import json, zlib, os, base64

""" Load variables from local file system. """
load_dotenv()

""" Set session configuration parameters. """
app.config["SESSION_TYPE"] = "filesystem"
app.config["SESSION_FILE_DIR"] = "./ptmvis/session/"
app.config["SESSION_KEY_PREFIX"] = "ptmvis:"
app.config["SESSION_USE_SIGNER"] = False
app.config["SESSION_PERMANENT"] = False
app.config["SESSION_FILE_THRESHOLD"] = 1000
app.config["MAX_CONTENT_LENGTH"] = 200 * 1024 * 1024  # Limit content lengths to 200 MB.

""" Set API parameters. """
api_parameters = {"URL": os.getenv("URL"), "DEBUG": os.getenv("DEBUG")}

""" Start session """
Session(app)

""" Definition of session keys """
MODIFICATIONS_DATA = "modifications_data_frame"
AVAILABLE_PROTEINS = "available_proteins"


def request_to_json(request_content: any) -> object:
    """Decompress zlib encoded JSON request into Python object.

    Uses UTF-8 decoding.

    Parameters
    ----------
    request_content : ReadableBuffer
        zlib encoded JSON request content.

    Returns
    -------
    object
        Python object yielding the content of the JSON request.
    """
    inflated_request_data = zlib.decompress(request_content)
    json_string_request_data = inflated_request_data.decode("utf8")
    return json.loads(json_string_request_data)


@app.route("/example_session", methods=["GET"])
def example_session():
    session.clear( )
    with open( "./ptmvis/static/resources/example_modifications_data.json", "r" ) as example_modifications_data :
        session[MODIFICATIONS_DATA] = json.load( example_modifications_data )
    return "Ok"

@app.route("/download_session", methods=["GET"])
def download_session():
    zlib_compress = zlib.compressobj( 6, zlib.DEFLATED, zlib.MAX_WBITS )
    compressed_session_bytes = zlib_compress.compress( bytes(json.dumps(session[MODIFICATIONS_DATA]), "utf-8") ) + zlib_compress.flush( )
    encoded_session = base64.b64encode( compressed_session_bytes ).decode("ascii")
    print( encoded_session )
    return encoded_session

@app.route("/restart_session", methods=["POST"])
def restart_session():
    session_data = zlib.decompress( base64.b64decode( request.data ) ).decode( )
    session[MODIFICATIONS_DATA] = json.loads(session_data)
    return "Ok"

@app.route("/process_search_engine_output", methods=["POST"])
def process_search_engine_output():
    """Route to process CSV data or SAGE, ionbot or MSFragger search engine output."""
    response = "Ok"
    try:
        json_request_data = request_to_json(request.data)
        content_type = json_request_data["contentType"]
        json_user_data = utils.parse_user_input(
            StringIO(json_request_data["content"]), content_type
        )
        session[MODIFICATIONS_DATA] = json_user_data
        if api_parameters["DEBUG"] :
            with open( "./dump.json", "w+" ) as dumpfile :
                dumpfile.write( json.dumps( session[MODIFICATIONS_DATA], indent = 3 ) )
    except TypeError as e:
        response = "Failed: " + str(e)
    finally:
        return response


@app.route("/get_available_proteins", methods=["GET"])
def get_available_proteins():
    """Route to retrieve all available proteins of the session as a JSON."""
    protein_entries = []
    if MODIFICATIONS_DATA in session:
        protein_identifiers = [_ for _ in session[MODIFICATIONS_DATA]["proteins"]]
        with open(
            "./ptmvis/static/resources/uniprot_id_mapping.json", "r"
        ) as protein_name_mapping_in:
            protein_name_mapping = json.load(protein_name_mapping_in)
        # Construct entry per protein in input data.
        for protein_identifier in protein_identifiers:
            if protein_identifier in protein_name_mapping:
                protein_name = protein_name_mapping[protein_identifier]
            else:
                protein_name = "N/A"
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
                    modifications.append(modification["modification_unimod_name"])
            modified_positions = sorted(modified_positions)
            modifications = list(set(modifications))
            protein_entry = {
                "id": protein_identifier,
                "name": protein_name,
                "modified_positions": len(modified_positions),
                "unique_modifications": len(modifications),
                "modifications": "$".join(modifications),
            }
            protein_entries.append(protein_entry)
    return protein_entries


@app.route("/get_overview_data", methods=["GET"])
def get_overview_data():
    """Route to retrieve nodes and links data for modifications overview graph, representing all modifications of a dataset."""
    modifications =  { } # Stores all present modifications together with their count.
    modification_co_occurrence = { } # Maps pairs of modification display names to their co-occurrence count.
    modification_classification_counts = { }
    if MODIFICATIONS_DATA in session:
        protein_identifiers = [_ for _ in session[MODIFICATIONS_DATA]["proteins"]]
        for protein_identifier in protein_identifiers:
            modification_data = session[MODIFICATIONS_DATA]["proteins"][
                protein_identifier
            ]["positions"]
            for position in modification_data:
                modifications_at_position = []
                for modification in modification_data[position][
                    "modifications"
                ]:
                    modification_unimod_name = modification["modification_unimod_name"]
                    modifications.setdefault( modification_unimod_name, modification )

                    modifications[ modification_unimod_name ].setdefault( "count", 0 )
                    modifications[ modification_unimod_name ][ "count" ] += 1

                    modifications[ modification_unimod_name ].setdefault( "occurrence", [ ] )
                    modifications[ modification_unimod_name ][ "occurrence" ].append( protein_identifier )

                    modifications_at_position.append(modification_unimod_name)

                    modification_classification = modification[ "modification_classification" ]
                    modification_classification_counts.setdefault( modification_classification, 0 )
                    modification_classification_counts[ modification_classification ] += 1

                for modifications_pair_tuple in combinations(set(modifications_at_position), 2):
                    modifications_pair = "@".join( sorted( list( modifications_pair_tuple ) ) )
                    modification_co_occurrence.setdefault( modifications_pair, 0 )
                    modification_co_occurrence[ modifications_pair ] += 1

        for modification_unimod_name, modification in modifications.items( ) :
            modification[ "frequency" ] = round( len( modification[ "occurrence" ] ) / len( protein_identifiers ), 4 )

        modifications_by_mass_shift = sorted( list( modifications.keys( ) ), key = lambda k : -modifications[k][ "mass_shift" ] if modifications[k][ "mass_shift" ] != "null" else 0.0 )
        modifications_by_count = sorted( list( modifications.keys( ) ), key = lambda k : -modifications[k][ "count" ] )

        return [
            modifications,
            [ modifications_by_mass_shift, modifications_by_count ],
            modification_co_occurrence,
            modification_classification_counts
        ]


@app.route("/get_protein_data", methods=["POST"])
def get_protein_data():
    response = {"status": "ok", "content": None}
    json_request_data = request_to_json(request.data)
    # Try to fetch PDB format structure for UniProt identifier from AlphaFold database.
    if json_request_data["pdb_text"] != None:
        structure, pdb_text = utils.parse_structure(json_request_data["pdb_text"])
    else:
        if "structure" in session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_id"] ] :
            pdb_text = utils._brotli_decompress( session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_id"] ]["structure"] )
            structure = utils.parse_structure( pdb_text )
        else :
            structure, pdb_text = utils.get_structure(json_request_data["uniprot_id"])
            # Extract protein sequence from structure and store it in session data.
            session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_id"] ]["sequence"] = utils.get_sequence_from_structure(structure)
            session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_id"] ][ "structure" ] = utils._brotly_compress(pdb_text)
    if structure != None:
        # Extract annotations for protein from UniProt.
        if not "annotation" in session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_id"] ] :
            annotation = _map_uniprot_identifiers(
                [json_request_data["uniprot_id"]], "UniProtKB"
            )
            # TODO: Clean annotation.
            session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_id"] ]["annotation"] = { }
            session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_id"] ]["annotation"] = annotation[ "results" ][ 0 ][ "to" ]

        # Compute contacts from structure and store them in session data.
        if not "contacts" in session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_id"] ] :
            session[MODIFICATIONS_DATA]["meta_data"]["distance_cutoff"] = int(json_request_data["cutoff"])
            session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_id"] ][ "contacts" ] = { }
            distance_matrix = utils.get_distance_matrix(structure)
            contacts = utils.get_contacts(
                distance_matrix,
                float(json_request_data["cutoff"]),
            )
            for source_position, contacts_list in contacts.items( ) :
                session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_id"] ][ "contacts" ][ source_position ] = [ ]
                for contact_position in contacts_list :
                    session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_id"] ][ "contacts" ][ source_position ].append( ( contact_position, distance_matrix[ source_position, contact_position ] ) )
        # Construct response.
        response[ "content" ] = deepcopy( session[MODIFICATIONS_DATA]["proteins"][ json_request_data["uniprot_id"] ] )
        response[ "content" ][ "structure" ] = utils._brotli_decompress( response[ "content" ][ "structure" ] )
    else :
        response = {"status": "Failed: No protein structure available.", "content": None}
    return response


def _map_uniprot_identifiers(identifiers, target_db):
    job_id = submit_id_mapping(
        from_db="UniProtKB_AC-ID", to_db=target_db, ids=identifiers
    )
    if check_id_mapping_results_ready(job_id):
        link = get_id_mapping_results_link(job_id)
        results = get_id_mapping_results_search(link)
    return results
